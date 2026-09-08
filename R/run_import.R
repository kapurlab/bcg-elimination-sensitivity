# Infected animals imported into the herd each year, under the paper's regime
# and under 18-month immunity with annual whole-herd revaccination.
# Run from the project root: Rscript R/run_import.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/profile.R")
herds <- load_herds("data"); R0 <- median(herds$R0); N <- 100
K <- 40; SL <- 0.1
regimes <- list(
  "Lifelong immunity, calves at birth (paper)" = list(prof = make_profile(K, SL, 30, 1, Inf), interval = NA, absorb = TRUE),
  "18-month immunity, annual whole-herd revaccination" = list(prof = make_profile(K, SL, 30, 1.5, 1.55), interval = 1, absorb = FALSE))
imports <- c(0, 0.25, 0.5, 1, 2, 5)     # infected animals entering per 100 head per year

res <- bind_rows(lapply(names(regimes), function(rn) bind_rows(lapply(imports, function(m) {
  r <- regimes[[rn]]
  s <- simulate_profile(R0, N, BASE$e_s, BASE$e_i, r$prof, stage_len = SL, interval = r$interval,
                        imports = m, absorb = r$absorb, horizon = 150, dt = 0.25)
  pv <- s$prev[, 1]
  data.frame(regime = rn, imports = m,
             T_1pct = crossing_time(s$times, pv, 1e-2), T_0.1pct = crossing_time(s$times, pv, 1e-3),
             prev_50y = pv[which.min(abs(s$times - 50))], prev_150y = pv[length(pv)],
             traj = I(list(data.frame(t = s$times, prev = pv))))
}))))
write.csv(res %>% select(-traj), "output/tables/imports.csv", row.names = FALSE)

tr <- res %>% select(regime, imports, traj) %>% tidyr::unnest(traj) %>% filter(t <= 60) %>%
  mutate(imports = factor(imports), regime = factor(regime, levels = names(regimes)))
g20 <- ggplot(tr, aes(t, prev, colour = imports)) +
  geom_hline(yintercept = c(1e-2, 1e-3), linetype = 3, colour = "grey55") +
  geom_text(data = data.frame(t = 60, prev = c(1e-2, 1e-3), lab = c("1%", "0.1%")), aes(t, prev, label = lab),
            inherit.aes = FALSE, hjust = 1, vjust = -0.3, size = 3, colour = "grey40") +
  geom_line(linewidth = 1) +
  scale_y_log10(breaks = c(0.5, 0.1, 0.01, 0.001, 1e-4), labels = c("50%", "10%", "1%", "0.1%", "0.01%"), limits = c(1e-4, 0.8)) +
  scale_colour_manual(values = c("#1F618D", "#7FB3D5", "#82E0AA", "#F9E79F", "#F1948A", "#7B1E3A"),
                      name = "Infected animals bought per 100 head per year") +
  facet_wrap(~regime) +
  labs(x = "Years since the start of vaccination", y = "Herd prevalence (log scale)",
       title = "Imports of infected animals set a floor under prevalence, median herd",
       subtitle = sprintf("R0 = %.1f, direct efficacy 58%%, indirect 74%%. Herd size held constant: each import replaces a birth.", R0)) +
  guides(colour = guide_legend(nrow = 1)) + theme_pub(12) + theme(plot.title.position = "plot")
ggsave("output/figs/fig20_imports.png", g20, width = 12, height = 6, dpi = 200, bg = "white")
cat("Imports done\n")
print(res %>% select(-traj) %>% mutate(across(c(T_1pct, T_0.1pct), ~ round(.x)), across(c(prev_50y, prev_150y), ~ signif(.x, 2))))
