# Infected purchases in a low-R0 setting with test-and-removal, as in the UK:
# within-herd R0 near 1.1 once annual testing removes most infected animals.
# Run from the project root: Rscript R/run_import_uk.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/profile.R")
N <- 100; K <- 40; SL <- 0.1
settings <- list(
  "Ethiopia: R0 2.8, no testing"                        = list(R0 = 2.808, v = 0),
  "UK-like: R0 1.1 with annual test-and-removal"        = list(R0 = 1.1, v = 0.7))
regimes <- list(
  "No vaccination"                                     = list(prof = make_profile(K, SL, 30, 1, Inf), interval = NA, absorb = TRUE, p = 0),
  "Lifelong immunity, calves at birth"                 = list(prof = make_profile(K, SL, 30, 1, Inf), interval = NA, absorb = TRUE, p = 1),
  "18-month immunity, annual whole-herd revaccination" = list(prof = make_profile(K, SL, 30, 1.5, 1.55), interval = 1, absorb = FALSE, p = 1))
imports <- c(0, 0.25, 0.5, 1, 2, 5)

res <- bind_rows(lapply(names(settings), function(sn) bind_rows(lapply(names(regimes), function(rn) bind_rows(lapply(imports, function(m) {
  st <- settings[[sn]]; r <- regimes[[rn]]
  s <- simulate_profile(st$R0, N, BASE$e_s, BASE$e_i, r$prof, stage_len = SL, interval = r$interval, p = r$p,
                        imports = m, absorb = r$absorb, v = st$v, horizon = 100, dt = 0.25)
  pv <- s$prev[, 1]
  data.frame(setting = sn, regime = rn, imports = m,
             prev_start = pv[1], prev_20y = pv[which.min(abs(s$times - 20))], prev_100y = pv[length(pv)],
             T_0.1pct = crossing_time(s$times, pv, 1e-3),
             traj = I(list(data.frame(t = s$times, prev = pv))))
}))))))
write.csv(res %>% select(-traj), "output/tables/imports_uk.csv", row.names = FALSE)

tr <- res %>% select(setting, regime, imports, traj) %>% unnest(traj) %>% filter(t <= 40) %>%
  mutate(imports = factor(imports), setting = factor(setting, levels = names(settings)), regime = factor(regime, levels = names(regimes)))
g21 <- ggplot(tr, aes(t, pmax(prev, 1e-5), colour = imports)) +
  geom_hline(yintercept = c(1e-2, 1e-3), linetype = 3, colour = "grey55") +
  geom_line(linewidth = 1) +
  scale_y_log10(breaks = c(0.5, 0.1, 0.01, 0.001, 1e-4), labels = c("50%", "10%", "1%", "0.1%", "0.01%"), limits = c(1e-4, 0.8)) +
  scale_colour_manual(values = c("#1F618D", "#7FB3D5", "#82E0AA", "#F9E79F", "#F1948A", "#7B1E3A"), name = "Infected animals bought per 100 head per year") +
  facet_grid(setting ~ regime, labeller = label_wrap_gen(30)) +
  labs(x = "Years since the start of the programme", y = "Herd prevalence (log scale)",
       title = "Infected purchases in a high-R0 setting without testing and a low-R0 setting with annual test-and-removal",
       subtitle = "UK-like herd: infected animals removed at 0.7 per year by testing, within-herd R0 = 1.1 with that removal included. Direct efficacy 58%, indirect 74%.") +
  guides(colour = guide_legend(nrow = 1)) + theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig21_imports_uk.png", g21, width = 13, height = 7.5, dpi = 200, bg = "white")
cat("UK imports done\n")
print(res %>% select(-traj) %>% mutate(across(c(prev_start, prev_20y, prev_100y), ~ signif(.x, 2)), T_0.1pct = round(T_0.1pct)) %>% as.data.frame())
