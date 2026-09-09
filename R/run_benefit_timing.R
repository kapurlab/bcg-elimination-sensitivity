# When does the benefit of vaccination arrive, relative to when herds become
# free of infection? A grid of probabilities cannot answer this because it has
# no time axis. Here the same strata are followed through time.
#
# For each herd size and R0, the burden remaining under vaccination is
# expressed as a share of the burden under no vaccination at the same moment.
# That falls from 1 towards a plateau. If the curves for different strata lie
# on top of one another, the timing of benefit does not depend on herd
# demography, which is the claim being tested.
#
# Run from the project root: Rscript R/run_benefit_timing.R (about 5 min)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/stochastic.R")
set.seed(20260909)
N_grid <- c(10, 20, 50, 100, 200); R0_grid <- c(1.5, 2, 3, 5)
REPS <- 2000; YEARS <- 50; TSTEP <- 60

grid <- expand.grid(N = N_grid, R0 = R0_grid)
CACHE <- "output/tables/benefit_timing.csv"
if (file.exists(CACHE)) {
  tr <- read.csv(CACHE); cat("using cached trajectories from", CACHE, "\n")
} else {
  tr <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    sq <- stoch_scenario(g$R0, g$N, p = 0, reps = REPS, years = YEARS, tstep = TSTEP, return_traj = TRUE)
    vx <- stoch_scenario(g$R0, g$N, D = 1.5, p = 1, reps = REPS, years = YEARS, tstep = TSTEP, return_traj = TRUE)
    data.frame(N = g$N, R0 = g$R0, yr = sq$yr,
               prev_sq = sq$mean_prev, prev_vx = vx$mean_prev, p_free_vx = vx$p_free)
  })) %>% mutate(remaining = prev_vx / prev_sq)
  write.csv(tr, CACHE, row.names = FALSE)
}

# summary: when is the burden reduction mostly achieved, and when are herds free
sm <- tr %>% group_by(N, R0) %>%
  summarise(plateau = min(remaining),
            yr_half_of_fall = yr[which(remaining <= 1 - 0.5 * (1 - min(remaining)))][1],
            yr_90_of_fall = yr[which(remaining <= 1 - 0.9 * (1 - min(remaining)))][1],
            yr_half_free = { i <- which(p_free_vx >= 0.5); if (length(i)) yr[i[1]] else NA },
            .groups = "drop")
write.csv(sm, "output/tables/benefit_timing_summary.csv", row.names = FALSE)

half_med <- median(sm$yr_half_of_fall, na.rm = TRUE)
never <- sum(is.na(sm$yr_half_free))
long <- tr %>%
  transmute(N, R0, yr,
            `Burden: infection remaining, relative to the same herd unvaccinated` = remaining,
            `Freedom: share of herds with no infected animal` = p_free_vx) %>%
  pivot_longer(-c(N, R0, yr), names_to = "metric", values_to = "value")
med <- long %>% group_by(metric, yr) %>% summarise(m = median(value), .groups = "drop")

g37 <- ggplot(long, aes(yr, value, group = interaction(N, R0))) +
  geom_vline(xintercept = half_med, linetype = 2, colour = "grey45") +
  geom_line(colour = "#7FB3D5", alpha = 0.5, linewidth = 0.55) +
  geom_line(data = med, aes(yr, m), inherit.aes = FALSE, colour = "#1B4F72", linewidth = 1.3) +
  facet_wrap(~metric, labeller = label_wrap_gen(46)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.02)) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  labs(x = "Years since the start of vaccination", y = NULL,
       title = "Benefit arrives early; freedom from infection arrives late or not at all",
       subtitle = sprintf(paste0("Pale lines: 20 strata, herd size 10 to 200 crossed with R0 1.5 to 5, 2000 herds each. Dark line: their median.\n",
                                 "An 18-month vaccine given to calves at birth, no boosters. Dashed line at year %.0f, when half of the\n",
                                 "total burden reduction has been achieved. In %d of the 20 strata, fewer than half of herds are ever free."),
                          half_med, never)) +
  theme_pub(12) + theme(plot.title.position = "plot")
ggsave("output/figs/fig37_benefit_timing.png", g37, width = 11, height = 6.2, dpi = 200, bg = "white")

cat("Benefit timing done\n")
print(as.data.frame(sm))
cat("\nmedian across strata: half the fall by year", half_med,
    "| 90% of the fall by year", median(sm$yr_90_of_fall, na.rm = TRUE),
    "| half of herds free by year", median(sm$yr_half_free, na.rm = TRUE),
    "(in", sum(is.na(sm$yr_half_free)), "of 20 strata that is never reached)\n")

# Figure 38 (the same strata drawn separately) is built by R/run_fig38.R,
# which reads the tables written above.
