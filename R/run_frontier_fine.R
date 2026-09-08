# Finer grid for the feasibility frontier, with the Ethiopian survey herds
# overlaid at their posterior-median R0 and tested herd size.
# Run from the project root: Rscript R/run_frontier_fine.R (about 15 min)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/revaccination.R")
set.seed(20260909)
herds <- load_herds("data")
R0_grid <- c(1.2, 1.35, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6.5, 8)
N_grid <- c(5, 8, 12, 20, 35, 50, 75, 100, 150, 200)
schedules <- list(
  "No vaccination (status quo)" = list(D = Inf, interval = NA, revacc = FALSE, p = 0),
  "Lifelong immunity, calves only (paper assumption)" = list(D = Inf, interval = NA, revacc = FALSE, p = 1),
  "18-month immunity, calves only" = list(D = 1.5, interval = NA, revacc = FALSE, p = 1),
  "18-month immunity, annual whole-herd booster" = list(D = 1.5, interval = 1, revacc = TRUE, p = 1))
REPS <- 200; YEARS <- 50
design <- expand.grid(R0 = R0_grid, N = N_grid, schedule = names(schedules), stringsAsFactors = FALSE)
res <- bind_rows(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, ]; sc <- schedules[[d$schedule]]
  et <- stoch_revacc(d$R0, d$N, BASE$e_s, BASE$e_i, D = sc$D, k = 20, p = sc$p,
                     interval = if (is.na(sc$interval)) 1 else sc$interval, reps = REPS, years = YEARS,
                     revacc = sc$revacc, redose_V = TRUE)
  data.frame(d, p_free_20y = mean(et <= 20), p_free_50y = mean(et <= 50))
}))
write.csv(res, "output/tables/frontier_fine.csv", row.names = FALSE)

prog_pal <- c("No vaccination (status quo)" = "grey45", "Lifelong immunity, calves only (paper assumption)" = "#7FB3D5",
              "18-month immunity, calves only" = "#F1948A", "18-month immunity, annual whole-herd booster" = "#1E8449")
gb <- res %>% pivot_longer(c(p_free_20y, p_free_50y), names_to = "horizon", values_to = "p") %>%
  mutate(horizon = recode(horizon, p_free_20y = "Chance of freedom by year 20", p_free_50y = "Chance of freedom by year 50"),
         lR = log10(R0), lN = log10(N), schedule = factor(schedule, levels = names(schedules)))
hp <- herds %>% mutate(lR = log10(pmin(R0, 8)), lN = log10(pmin(herdsize, 200)))
g31 <- ggplot(gb, aes(lR, lN)) +
  geom_point(data = hp, aes(lR, lN), inherit.aes = FALSE, colour = "grey30", fill = "#F9E79F", shape = 21, size = 2.2, alpha = 0.9) +
  geom_contour(aes(z = p, colour = schedule), breaks = 0.5, linewidth = 1.1) +
  geom_contour(aes(z = p, colour = schedule), breaks = 0.8, linewidth = 1.1, linetype = 2) +
  scale_colour_manual(values = prog_pal, name = NULL, labels = function(x) stringr::str_wrap(x, 30)) +
  scale_x_continuous(breaks = log10(c(1.2, 1.5, 2, 3, 5, 8)), labels = c(1.2, 1.5, 2, 3, 5, 8)) +
  scale_y_continuous(breaks = log10(c(5, 10, 20, 50, 100, 200)), labels = c(5, 10, 20, 50, 100, 200)) +
  facet_wrap(~horizon) +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "Feasibility frontier with the 57 Ethiopian survey herds overlaid",
       subtitle = "Lines: herd size and R0 at which a programme gives a 50% (solid) or 80% (dashed) chance of no infected animal; 200 replicates per point, 11 by 10 grid.\nYellow points: survey herds at posterior-median R0 and tested size (values above 8 or 200 drawn at the edge). Herds below and left of a line are within reach.") +
  guides(colour = guide_legend(ncol = 4)) + theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig31_frontier_with_herds.png", g31, width = 14, height = 7.5, dpi = 200, bg = "white")
cat("Fine frontier done\n")
