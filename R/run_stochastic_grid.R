# Stochastic (SimInf) results on a conceptual grid of within-herd R0 and herd
# size, so that results are not tied to the Ethiopian survey herds. Uses the
# authors' transition list with k-stage near-fixed protection and whole-herd
# campaign pulses (stoch_revacc in R/revaccination.R). Elimination is the
# first day with no infected animal in the herd.
# Run from the project root: Rscript R/run_stochastic_grid.R (about 20 min)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/revaccination.R")
set.seed(20260908)
R0_grid <- c(1.2, 1.5, 2, 3, 5, 8)
N_grid <- c(5, 10, 20, 50, 100, 200)
schedules <- list(
  "No vaccination (status quo)"                       = list(D = Inf, interval = NA, revacc = FALSE, p = 0),
  "Lifelong immunity, calves only (paper assumption)" = list(D = Inf, interval = NA, revacc = FALSE),
  "18-month immunity, calves only"                   = list(D = 1.5, interval = NA, revacc = FALSE),
  "18-month immunity, annual whole-herd booster"     = list(D = 1.5, interval = 1, revacc = TRUE),
  "18-month immunity, six-monthly whole-herd booster" = list(D = 1.5, interval = 0.5, revacc = TRUE))
REPS <- 200; YEARS <- 100

design <- expand.grid(R0 = R0_grid, N = N_grid, schedule = names(schedules), stringsAsFactors = FALSE)
res <- bind_rows(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, ]; sc <- schedules[[d$schedule]]
  pv <- if (is.null(sc$p)) 1 else sc$p
  et <- stoch_revacc(d$R0, d$N, BASE$e_s, BASE$e_i, D = sc$D, k = 20, p = pv,
                     interval = if (is.na(sc$interval)) 1 else sc$interval, reps = REPS, years = YEARS,
                     revacc = sc$revacc, redose_V = TRUE)
  data.frame(d, p_free_10y = mean(et <= 10), p_free_20y = mean(et <= 20), p_free_50y = mean(et <= 50),
             p_free_100y = mean(is.finite(et)), median_T = median(et),
             Rv = R_vacc(d$R0, BASE$e_s, BASE$e_i, sc$D, pv))
}))
write.csv(res, "output/tables/stochastic_grid.csv", row.names = FALSE)

# --- figure 22: probability free by year 20 ------------------------------------
lab <- function(x) sprintf("%.0f%%", 100 * x)
d <- res %>% mutate(schedule = factor(schedule, levels = names(schedules)),
                    R0f = factor(R0), Nf = factor(N))
g22 <- ggplot(d, aes(R0f, Nf, fill = p_free_20y)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = lab(p_free_20y), colour = p_free_20y > 0.55), size = 3.4, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(0, 1), labels = scales::percent,
                      name = "Probability the herd has no infected animal by year 20") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  facet_wrap(~schedule, ncol = 3, labeller = label_wrap_gen(40)) +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "Stochastic elimination on a conceptual grid: herd size against R0, status quo and four programmes",
       subtitle = sprintf("%d replicate herds per cell, each starting at its endemic level; every calf dosed at birth. Direct efficacy 58%%, indirect 74%%, turnover 0.27 per year.\nEthiopian survey herds: R0 median 2.8, quartiles 2.0 and 4.1; sizes 21 to 146.", REPS)) +
  theme_pub(11) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig22_stochastic_grid_20y.png", g22, width = 15, height = 8, dpi = 200, bg = "white")

# --- figure 23: median years to freedom ---------------------------------------
d2 <- d %>% mutate(never = !is.finite(median_T), fill_val = ifelse(never, NA, pmin(median_T, 100)),
                   label = ifelse(never, "> 100", sprintf("%.0f", median_T)))
g23 <- ggplot(d2, aes(R0f, Nf, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label, colour = never | fill_val > 60), size = 3.4, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(0, 100), na.value = "grey88",
                      name = "Median years to no infected animal (grey: more than 100 in half of herds)") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  facet_wrap(~schedule, ncol = 3, labeller = label_wrap_gen(40)) +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "Median years to a herd with no infected animal",
       subtitle = "Read each programme against the status quo panel: small herds fade out by chance without any vaccine; large herds need the reproduction number under vaccination below 1.") +
  theme_pub(11) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig23_stochastic_grid_median.png", g23, width = 15, height = 8, dpi = 200, bg = "white")
cat("Stochastic grid done\n")
print(res %>% filter(grepl("annual", schedule)) %>% mutate(v = lab(p_free_20y)) %>% select(R0, N, v) %>% pivot_wider(names_from = R0, values_from = v))
