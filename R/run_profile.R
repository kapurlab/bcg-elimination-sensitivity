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
  "Plateau to 6 months, gone by 36"             = make_profile(K, SL, 30, 0.5, 3.0),
  "Rise to 2 months, plateau to 9, gone by 18"  = make_profile(K, SL, 60, 0.75, 1.5))
sched <- data.frame(schedule = c("No booster", "Booster every 6 months", "Booster annually", "Booster every 24 months"),
                    interval = c(NA, 0.5, 1, 2))

# protection-years per dose (area under the curve) for reference
area <- sapply(profiles, function(pr) sum(pr) * SL)

design <- merge(data.frame(profile = names(profiles), area = area), sched)
design$T <- mapply(function(pn, iv) {
  s <- simulate_profile(R0, N, BASE$e_s, BASE$e_i, prof_p = profiles[[pn]], stage_len = SL, interval = iv)
  crossing_time(s$times, s$prev[, 1], 1e-3)
}, design$profile, design$interval)

# uncertainty: joint posterior draws of direct and indirect efficacy from the
# DST1 fit; timing of the profile is held fixed (no data on it yet)
set.seed(20260908)
post <- load_efficacy_posterior("data")
n_draw <- 40
dr <- post[sample(nrow(post), n_draw), c("eff_S", "eff_I")]
es_q <- quantile(post$eff_S, c(0.025, 0.5, 0.975))
Tdraws <- mapply(function(pn, iv) {
  vapply(seq_len(n_draw), function(i) {
    s <- simulate_profile(R0, N, dr$eff_S[i], dr$eff_I[i], prof_p = profiles[[pn]], stage_len = SL, interval = iv)
    crossing_time(s$times, s$prev[, 1], 1e-3)
  }, numeric(1))
}, design$profile, design$interval, SIMPLIFY = FALSE)
design$T_med <- sapply(Tdraws, median)
design$T_lo <- sapply(Tdraws, quantile, 0.025)
design$T_hi <- sapply(Tdraws, quantile, 0.975)
design$share_never <- sapply(Tdraws, function(x) mean(!is.finite(x)))
write.csv(design, "output/tables/profile_schedule.csv", row.names = FALSE)
base_T <- T_single(R0, BASE$e_s, BASE$e_i, Inf, 1, N = N)

# --- figure 17: the profiles with a credible band from the efficacy posterior ---
pp <- bind_rows(lapply(names(profiles), function(n) data.frame(profile = n, t = tmid * 12, rel = profiles[[n]])))
pp$profile <- factor(pp$profile, levels = names(profiles))
pp <- pp %>% mutate(mid = rel * es_q[2], lo = rel * es_q[1], hi = rel * es_q[3])
g17 <- ggplot(pp, aes(t, mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = profile), alpha = 0.25) +
  geom_line(aes(colour = profile), linewidth = 1) +
  facet_wrap(~profile, ncol = 3, labeller = label_wrap_gen(34)) +
  scale_colour_manual(values = c("#F1948A", "#7FB3D5", "#1F618D", "#82E0AA", "#BB8FCE", "#F5B041"), guide = "none") +
  scale_fill_manual(values = c("#F1948A", "#7FB3D5", "#1F618D", "#82E0AA", "#BB8FCE", "#F5B041"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 48, 12)) + scale_y_continuous(labels = scales::percent, limits = c(0, 0.75)) +
  labs(x = "Months since dose", y = "Reduction in susceptibility",
       title = "Protection profiles with the 95% credible band of the direct-efficacy plateau",
       subtitle = sprintf("Plateau at the DST1 posterior median (%.0f%%), band %.0f%% to %.0f%%. Timing of rise, plateau and decline is assumed, not estimated.\nProtection-years per dose at the median, in panel order: %s",
                          100 * es_q[2], 100 * es_q[1], 100 * es_q[3], paste(sprintf("%.2f", area), collapse = ", "))) +
  theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig17_protection_profiles.png", g17, width = 12, height = 7, dpi = 200, bg = "white")

# --- figure 18: tiles with credible intervals ------------------------------
fmt <- function(x) ifelse(is.finite(x), sprintf("%.0f", x), "never")
d <- design %>% mutate(profile = factor(profile, levels = names(profiles)),
                       schedule = factor(schedule, levels = rev(sched$schedule)),
                       never = !is.finite(T_med), fill_val = ifelse(never, NA, pmin(T_med, 120)),
                       label = ifelse(never, sprintf("never\n(%.0f%% of draws)", 100 * share_never),
                                      sprintf("%s\n(%s to %s; %.0f%% never)", fmt(T_med), fmt(T_lo), fmt(T_hi), 100 * share_never)))
g18 <- ggplot(d, aes(profile, schedule, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = never | fill_val > 80), size = 3.3, fontface = "bold", lineheight = 0.9) +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(25, 120), na.value = "grey88",
                      name = "Years to herd prevalence < 0.1% (posterior median)", breaks = c(25, 50, 75, 100, 120), labels = c("25", "50", "75", "100", "120+")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, 16)) +
  labs(x = "Protection profile after each dose", y = "Revaccination interval",
       title = "Schedule against protection profile, median herd, with 95% credible intervals",
       subtitle = sprintf("Intervals from %d joint posterior draws of direct and indirect efficacy (DST1 fit); profile timing fixed. Every calf dosed at birth; boosters to the whole herd. Lifelong case: %.0f years.", n_draw, base_T)) +
  theme_pub(12) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig18_profile_tiles.png", g18, width = 13, height = 6.5, dpi = 200, bg = "white")
cat("Profile done\n"); print(design %>% mutate(cell = sprintf("%s (%s-%s)", fmt(T_med), fmt(T_lo), fmt(T_hi))) %>% select(profile, schedule, cell) %>% pivot_wider(names_from = schedule, values_from = cell), width = 250)
