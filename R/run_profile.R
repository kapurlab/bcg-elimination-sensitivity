# Protection that rises, plateaus and declines with time since dose, against
# revaccination schedule. Run from the project root: Rscript R/run_profile.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/profile.R")
herds <- load_herds("data"); R0 <- median(herds$R0); N <- 100
K <- 40; SL <- 0.1                     # 40 stages of 0.1 year: the curve spans 4 years after a dose
tmid <- (seq_len(K) - 0.5) * SL

profiles <- list(
  "Step: full to 18 months, then none"          = make_profile(K, SL, 30, 1.5, 1.55),
  "Plateau to 12 months, gone by 24"            = make_profile(K, SL, 30, 1.0, 2.0),
  "Plateau to 12 months, gone by 36"            = make_profile(K, SL, 30, 1.0, 3.0),
  "Plateau to 18 months, gone by 36"            = make_profile(K, SL, 30, 1.5, 3.0),
  "Plateau to 6 months, gone by 36"             = make_profile(K, SL, 30, 0.5, 3.0))
sched <- data.frame(schedule = c("No booster", "Booster every 6 months", "Booster annually", "Booster every 24 months"),
                    interval = c(NA, 0.5, 1, 2))

# protection-years per dose (area under the curve) for reference
area <- sapply(profiles, function(pr) sum(pr) * SL)

design <- merge(data.frame(profile = names(profiles), area = area), sched)
design$T <- mapply(function(pn, iv) {
  s <- simulate_profile(R0, N, BASE$e_s, BASE$e_i, prof_p = profiles[[pn]], stage_len = SL, interval = iv)
  crossing_time(s$times, s$prev[, 1], 1e-3)
}, design$profile, design$interval)
write.csv(design, "output/tables/profile_schedule.csv", row.names = FALSE)
base_T <- T_single(R0, BASE$e_s, BASE$e_i, Inf, 1, N = N)

# --- figure 17: the profiles ------------------------------------------------
pp <- bind_rows(lapply(names(profiles), function(n) data.frame(profile = n, t = tmid * 12, prot = profiles[[n]] * BASE$e_s)))
pp$profile <- factor(pp$profile, levels = names(profiles))
g17 <- ggplot(pp, aes(t, prot, colour = profile)) +
  geom_step(linewidth = 1) +
  scale_colour_manual(values = c("#F1948A", "#7FB3D5", "#1F618D", "#82E0AA", "#BB8FCE"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 48, 6)) + scale_y_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  labs(x = "Months since dose", y = "Reduction in susceptibility",
       title = "Protection profiles compared: a 30-day rise, a plateau at the paper's 58%, then decline",
       subtitle = sprintf("Protection-years per dose (area under each curve), left to right in the legend order: %s", paste(sprintf("%.2f", area), collapse = ", "))) +
  guides(colour = guide_legend(nrow = 2)) + theme_pub(12) + theme(plot.title.position = "plot")
ggsave("output/figs/fig17_protection_profiles.png", g17, width = 11, height = 5.5, dpi = 200, bg = "white")

# --- figure 18: tiles -----------------------------------------------------------
d <- design %>% mutate(profile = factor(profile, levels = names(profiles)),
                       schedule = factor(schedule, levels = rev(sched$schedule)),
                       never = !is.finite(T), fill_val = ifelse(never, NA, pmin(T, 120)),
                       label = ifelse(never, "never", sprintf("%.0f", T)))
g18 <- ggplot(d, aes(profile, schedule, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = never | fill_val > 80), size = 4.6, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(25, 120), na.value = "grey88",
                      name = "Years to herd prevalence < 0.1%", breaks = c(25, 50, 75, 100, 120), labels = c("25", "50", "75", "100", "120+")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, 18)) +
  labs(x = "Protection profile after each dose", y = "Revaccination interval",
       title = "Gradual decline softens the diagonal: schedule against protection profile, median herd",
       subtitle = sprintf("Every calf dosed at birth; boosters to the whole herd, each dose restarting the profile. Paper's lifelong case: %.0f years. R0 = %.1f.", base_T, R0)) +
  theme_pub(12) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig18_profile_tiles.png", g18, width = 12, height = 6, dpi = 200, bg = "white")
cat("Profile done\n"); print(design %>% mutate(T = round(T)) %>% select(profile, schedule, T) %>% pivot_wider(names_from = schedule, values_from = T))
