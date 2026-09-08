# Alternative views of the stochastic results: difference tiles, feasibility
# boundary map, archetype slices, survival curves, and prevalence bands for
# the replacement policy. Run from the project root: Rscript R/run_alt_views.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/revaccination.R"); source("R/herd.R")
set.seed(20260908)
grid <- read.csv("output/tables/stochastic_grid.csv")
sched_levels <- unique(grid$schedule)
prog_pal <- c("No vaccination (status quo)" = "grey45",
              "Lifelong immunity, calves only (paper assumption)" = "#7FB3D5",
              "18-month immunity, calves only" = "#F1948A",
              "18-month immunity, annual whole-herd booster" = "#1E8449",
              "18-month immunity, six-monthly whole-herd booster" = "#82E0AA")

# --- fig 26: difference tiles, programme minus status quo ------------------------
ctrl <- grid %>% filter(schedule == sched_levels[1]) %>% select(R0, N, c20 = p_free_20y, c50 = p_free_50y)
diff <- grid %>% filter(schedule != sched_levels[1]) %>% left_join(ctrl, by = c("R0", "N")) %>%
  mutate(d20 = p_free_20y - c20, d50 = p_free_50y - c50, schedule = factor(schedule, levels = sched_levels[-1]))
g26 <- ggplot(diff, aes(factor(R0), factor(N), fill = d20)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.0f", 100 * d20), colour = abs(d20) > 0.5), size = 3.3, fontface = "bold") +
  scale_fill_gradient2(low = "#B7950B", mid = "white", high = "#7B1E3A", midpoint = 0, limits = c(-1, 1),
                       labels = function(x) sprintf("%+.0f", 100 * x), name = "Percentage points added to the probability of freedom by year 20") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  facet_wrap(~schedule, ncol = 2, labeller = label_wrap_gen(45)) +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "What vaccination adds: probability of freedom by year 20, programme minus status quo",
       subtitle = "Positive values are the vaccine's contribution. Small low-R0 herds gain little because they clear infection by chance anyway.") +
  theme_pub(11) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.8, "cm"))
ggsave("output/figs/fig26_difference_tiles.png", g26, width = 12, height = 9, dpi = 200, bg = "white")

# --- fig 27: feasibility boundary map ------------------------------------------------
gb <- grid %>% mutate(lR = log10(R0), lN = log10(N), schedule = factor(schedule, levels = sched_levels))
g27 <- ggplot(gb, aes(lR, lN, z = p_free_20y, colour = schedule)) +
  geom_contour(breaks = 0.5, linewidth = 1.1) +
  geom_contour(breaks = 0.8, linewidth = 1.1, linetype = 2) +
  scale_colour_manual(values = prog_pal, name = NULL, labels = function(x) stringr::str_wrap(x, 32)) +
  scale_x_continuous(breaks = log10(c(1.2, 1.5, 2, 3, 5, 8)), labels = c(1.2, 1.5, 2, 3, 5, 8)) +
  scale_y_continuous(breaks = log10(c(5, 10, 20, 50, 100, 200)), labels = c(5, 10, 20, 50, 100, 200)) +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "Feasibility frontier: where each programme gives a 50% (solid) or 80% (dashed) chance of freedom by year 20",
       subtitle = "Below and left of a line the programme reaches that chance; the frontier moves up and right as the programme improves.\nContours are interpolated from the 6 by 6 grid.") +
  guides(colour = guide_legend(ncol = 3)) + theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig27_feasibility_frontier.png", g27, width = 11, height = 8, dpi = 200, bg = "white")

# --- fig 28: archetype slices ----------------------------------------------------------
arch <- tribble(~archetype, ~N, ~R0,
  "Smallholder, 5 head, R0 2", 5, 2, "Smallholder, 10 head, R0 3", 10, 3,
  "Mid-size, 20 head, R0 3", 20, 3, "Commercial, 100 head, R0 2", 100, 2, "Commercial, 100 head, R0 3", 100, 3)
sl <- arch %>% left_join(grid, by = c("N", "R0")) %>%
  pivot_longer(c(p_free_20y, p_free_50y), names_to = "horizon", values_to = "p") %>%
  mutate(horizon = recode(horizon, p_free_20y = "by year 20", p_free_50y = "by year 50"),
         archetype = factor(archetype, levels = arch$archetype), schedule = factor(schedule, levels = sched_levels))
g28 <- ggplot(sl, aes(p, archetype, colour = schedule)) +
  geom_line(aes(group = archetype), colour = "grey80", linewidth = 2) +
  geom_point(size = 3.6, position = position_dodge(width = 0)) +
  scale_colour_manual(values = prog_pal, name = NULL, labels = function(x) stringr::str_wrap(x, 30)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~horizon) +
  labs(x = "Probability the herd has no infected animal", y = NULL,
       title = "Archetype slices: five representative herds, five programmes, two horizons",
       subtitle = "Grey bar spans the programmes for each herd; the status quo dot is the control.") +
  guides(colour = guide_legend(ncol = 3)) + theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig28_archetype_slices.png", g28, width = 12, height = 6, dpi = 200, bg = "white")

# --- fig 29: survival curves for archetype cells (re-simulated) -----------------------
schedules <- list(
  "No vaccination (status quo)" = list(D = Inf, interval = NA, revacc = FALSE, p = 0),
  "Lifelong immunity, calves only (paper assumption)" = list(D = Inf, interval = NA, revacc = FALSE, p = 1),
  "18-month immunity, calves only" = list(D = 1.5, interval = NA, revacc = FALSE, p = 1),
  "18-month immunity, annual whole-herd booster" = list(D = 1.5, interval = 1, revacc = TRUE, p = 1),
  "18-month immunity, six-monthly whole-herd booster" = list(D = 1.5, interval = 0.5, revacc = TRUE, p = 1))
cells <- tribble(~cell, ~N, ~R0, "5 head, R0 2", 5, 2, "10 head, R0 3", 10, 3, "20 head, R0 3", 20, 3, "100 head, R0 2", 100, 2, "100 head, R0 3", 100, 3, "200 head, R0 3", 200, 3)
yrs <- seq(0, 60, by = 1)
surv <- bind_rows(lapply(seq_len(nrow(cells)), function(i) bind_rows(lapply(names(schedules), function(sn) {
  sc <- schedules[[sn]]
  et <- stoch_revacc(cells$R0[i], cells$N[i], BASE$e_s, BASE$e_i, D = sc$D, k = 20, p = sc$p,
                     interval = if (is.na(sc$interval)) 1 else sc$interval, reps = 300, years = 60,
                     revacc = sc$revacc, redose_V = TRUE)
  data.frame(cell = cells$cell[i], schedule = sn, yr = yrs, p_free = sapply(yrs, function(t) mean(et <= t)))
}))))
write.csv(surv, "output/tables/survival_curves.csv", row.names = FALSE)
surv <- surv %>% mutate(cell = factor(cell, levels = cells$cell), schedule = factor(schedule, levels = sched_levels))
g29 <- ggplot(surv, aes(yr, p_free, colour = schedule)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = prog_pal, name = NULL, labels = function(x) stringr::str_wrap(x, 30)) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~cell, ncol = 3) +
  labs(x = "Years since the start of the programme", y = "Probability the herd has no infected animal",
       title = "Time course of freedom from infection in six representative herds, 300 replicates each",
       subtitle = "'Free by year 20' in the tile figures is one vertical slice through these curves.") +
  guides(colour = guide_legend(ncol = 3)) + theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig29_survival_curves.png", g29, width = 12, height = 8, dpi = 200, bg = "white")

# --- fig 30: prevalence bands under replacement policies (re-simulated) --------------
mod <- build_herd_model(k = 10)
pol <- tribble(~policy, ~h, ~pi_src,
  "All replacements home-bred", 1, 0, "Half purchased, source 5%", 0.5, 0.05,
  "Half purchased, source 20%", 0.5, 0.20, "All purchased, source 20%", 0, 0.20)
prog <- tribble(~programme, ~p, ~interval, "No vaccination", 0, NA, "18-month immunity, annual whole-herd booster", 1, 1)
hs <- crossing(NA_adults = c(5, 20, 100), R0 = 3, pol, prog)
traj <- bind_rows(lapply(seq_len(nrow(hs)), function(i) {
  d <- hs[i, ]
  tr <- run_herd(mod, R0 = d$R0, NA_adults = d$NA_adults, p = d$p, interval = d$interval, h = d$h, pi_src = d$pi_src,
                 reps = 200, years = 50, return_traj = TRUE)
  data.frame(d, tr)
}))
write.csv(traj, "output/tables/replacement_trajectories.csv", row.names = FALSE)
tp <- traj %>% mutate(policy = factor(policy, levels = pol$policy), programme = factor(programme, levels = prog$programme),
                      size = factor(sprintf("%d adults", NA_adults), levels = sprintf("%d adults", c(5, 20, 100))))
g30 <- ggplot(tp, aes(yr, prev_med, colour = programme, fill = programme)) +
  geom_ribbon(aes(ymin = prev_q25, ymax = prev_q75), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("grey40", "#1E8449"), name = NULL) + scale_fill_manual(values = c("grey40", "#1E8449"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(policy ~ size, labeller = label_wrap_gen(24)) +
  labs(x = "Years since the start of the programme", y = "Herd prevalence (median and interquartile band across 200 herds)",
       title = "Prevalence under each replacement policy, with and without vaccination, R0 = 3",
       subtitle = "Purchases from an infected source set the floor; vaccination lowers the level but cannot remove the floor.") +
  theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig30_replacement_prevalence_bands.png", g30, width = 12, height = 10, dpi = 200, bg = "white")
cat("Alternative views done\n")
