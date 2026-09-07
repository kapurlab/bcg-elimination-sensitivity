# Duration of protection under three vaccination strategies:
#   1. calves at birth only (the paper's design, as in run_oat.R)
#   2. calves at birth plus an annual campaign that vaccinates every
#      unprotected uninfected animal (S only)
#   3. calves at birth plus an annual campaign that vaccinates every uninfected
#      animal, restarting the protection clock of those still protected
# each under gradual (exponential, k = 1) and near-fixed (Erlang k = 20)
# waning. Run from the project root: Rscript R/run_revaccination.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R")
source("R/model.R")
source("R/revaccination.R")
set.seed(20260907)
herds <- load_herds("data")
qR0 <- quantile(herds$R0, c(0.25, 0.5, 0.75, 0.9))
rep_herds <- data.frame(label = sprintf(c("R0 25th pct (%.1f)", "R0 median (%.1f)",
                                          "R0 75th pct (%.1f)", "R0 90th pct (%.1f)"), qR0),
                        R0 = as.numeric(qR0))
VE <- 1 - (1 - BASE$e_s) * (1 - BASE$e_i)

D_grid <- c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 15, 20, Inf)
strategies <- data.frame(
  strategy = c("Calves at birth only", "Annual campaign: unprotected animals",
               "Annual campaign: all uninfected animals"),
  revacc = c(FALSE, TRUE, TRUE), redose_V = c(FALSE, FALSE, TRUE))
shapes <- data.frame(shape = c("Gradual waning (exponential)", "Near-fixed duration (Erlang k = 20)"),
                     k = c(1, 20))
design <- merge(merge(strategies, shapes), data.frame(D = D_grid))

# --- time-averaged vaccinated fraction and average R_v -----------------------
vf <- design %>% rowwise() %>%
  mutate(f_bar = vacc_fraction_revacc(D, k = k, revacc = revacc, redose_V = redose_V)) %>%
  ungroup()
vf_long <- vf %>% crossing(rep_herds) %>% mutate(Rv_bar = R0 * (1 - VE * f_bar))
write.csv(vf_long, "output/tables/revacc_vaccinated_fraction.csv", row.names = FALSE)

# --- deterministic time to elimination, representative herds ---------------
det <- bind_rows(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, ]
  sim <- simulate_revacc(rep_herds$R0, rep(100, 4), BASE$e_s, BASE$e_i, D = d$D, k = d$k,
                         revacc = d$revacc, redose_V = d$redose_V, horizon = 200, dt = 0.25)
  data.frame(d, herd = rep_herds$label, R0 = rep_herds$R0,
             T_elim = vapply(1:4, function(j) crossing_time(sim$times, sim$prev[, j], 1e-3), numeric(1)))
}))
write.csv(det, "output/tables/revacc_deterministic.csv", row.names = FALSE)

# --- deterministic, whole 57-herd population, 50-year horizon --------------
pop <- bind_rows(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, ]
  sim <- simulate_revacc(herds$R0, herds$herdsize, BASE$e_s, BASE$e_i, D = d$D, k = d$k,
                         revacc = d$revacc, redose_V = d$redose_V, horizon = 50, dt = 0.5)
  Th <- vapply(seq_len(nrow(herds)), function(j) crossing_time(sim$times, sim$prev[, j], 1e-3), numeric(1))
  data.frame(d, share_herds_free_50y = mean(Th <= 50),
             pop_prev_50y = sum(sim$infected[nrow(sim$infected), ]) / sum(herds$herdsize))
}))
write.csv(pop, "output/tables/revacc_population.csv", row.names = FALSE)

# --- stochastic check, median herd -----------------------------------------
D_stoch <- c(0.5, 1, 1.5, 2, 3, 5, 10, Inf)
sdesign <- merge(merge(strategies, shapes), data.frame(D = D_stoch))
sto <- bind_rows(lapply(seq_len(nrow(sdesign)), function(i) {
  d <- sdesign[i, ]
  et <- stoch_revacc(qR0[2], 44, BASE$e_s, BASE$e_i, D = d$D, k = d$k, reps = 300,
                     revacc = d$revacc, redose_V = d$redose_V)
  data.frame(d, median_T = median(et), q25 = quantile(et, 0.25), q75 = quantile(et, 0.75),
             share_eliminated_50y = mean(et <= 50), share_eliminated_200y = mean(is.finite(et)))
}))
write.csv(sto, "output/tables/revacc_stochastic_median_herd.csv", row.names = FALSE)

# --- figures -------------------------------------------------------------------
strat_pal <- c("Calves at birth only" = "#F1948A",
               "Annual campaign: unprotected animals" = "#7FB3D5",
               "Annual campaign: all uninfected animals" = "#82E0AA")
month_breaks <- c(3, 6, 12, 24, 60, 120, 240, 480)
to_months <- function(D) ifelse(is.infinite(D), 480, D * 12)   # 480 = lifelong on the axis
month_labels <- c("3", "6", "12", "24", "60", "120", "240", "life")

dd <- det %>% filter(herd %in% rep_herds$label[2:3]) %>%
  mutate(months = to_months(D), T_plot = ifelse(is.finite(T_elim), T_elim, NA),
         strategy = factor(strategy, levels = strategies$strategy))
g9 <- ggplot(dd, aes(months, T_plot, colour = strategy)) +
  geom_hline(yintercept = 50, linetype = 3, colour = "grey55") +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_x_log10(breaks = month_breaks, labels = month_labels) +
  scale_colour_manual(values = strat_pal, name = NULL) +
  facet_grid(herd ~ shape) +
  labs(x = "Duration of protection (months)", y = "Years to herd prevalence < 0.1%",
       title = "Duration of protection with and without an annual revaccination campaign",
       subtitle = "e_s 0.58, e_i 0.74, all calves vaccinated at birth. Gaps: elimination not reached within 200 years.\nUnder exponential waning re-dosing a protected animal changes nothing (memoryless), so the two campaign lines coincide.") +
  guides(colour = guide_legend(nrow = 1)) +
  theme_pub()
ggsave("output/figs/fig9_revaccination_duration.png", g9, width = 11, height = 7.5, dpi = 200, bg = "white")

need <- data.frame(herd = rep_herds$label, f_need = (1 - 1 / rep_herds$R0) / VE)
vv <- vf %>% mutate(months = to_months(D), strategy = factor(strategy, levels = strategies$strategy))
g10 <- ggplot(vv, aes(months, f_bar, colour = strategy)) +
  geom_hline(data = need, aes(yintercept = f_need), linetype = 2, colour = "grey45") +
  geom_text(data = need, aes(x = 3, y = f_need, label = herd), inherit.aes = FALSE,
            hjust = 0, vjust = -0.3, size = 3, colour = "grey30") +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_x_log10(breaks = month_breaks, labels = month_labels) +
  scale_colour_manual(values = strat_pal, name = NULL) +
  facet_wrap(~shape) +
  labs(x = "Duration of protection (months)", y = "Time-averaged protected fraction of the herd",
       title = "What the strategies deliver: the protected fraction of the herd",
       subtitle = "Dashed lines: fraction needed for R_v < 1 in herds at each R0 percentile, with 89% total efficacy.\nUnder exponential waning the two campaign lines coincide.") +
  guides(colour = guide_legend(nrow = 1)) +
  theme_pub()
ggsave("output/figs/fig10_revaccination_protected_fraction.png", g10, width = 11, height = 5.5, dpi = 200, bg = "white")

ss <- sto %>% mutate(months = to_months(D), strategy = factor(strategy, levels = strategies$strategy),
                     median_T = ifelse(is.finite(median_T), median_T, NA),
                     q75 = ifelse(is.finite(q75), q75, NA), q25 = ifelse(is.finite(q25), q25, NA))
g11 <- ggplot(ss, aes(months, median_T, colour = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_x_log10(breaks = month_breaks, labels = month_labels) +
  scale_colour_manual(values = strat_pal, name = NULL) + scale_fill_manual(values = strat_pal, name = NULL) +
  facet_wrap(~shape) +
  labs(x = "Duration of protection (months)", y = "Years to zero infected animals (median, IQR)",
       title = sprintf("Stochastic check, median herd (R0 = %.2f, N = 44), 300 replicates", qR0[2]),
       subtitle = "Gaps: fewer than half of replicates eliminated within 200 years") +
  guides(colour = guide_legend(nrow = 1), fill = "none") +
  theme_pub()
ggsave("output/figs/fig11_revaccination_stochastic.png", g11, width = 11, height = 5.5, dpi = 200, bg = "white")

cat("Revaccination done\n")
print(det %>% filter(herd == rep_herds$label[2]) %>% select(strategy, shape, D, T_elim) %>%
        pivot_wider(names_from = shape, values_from = T_elim), n = 60)
print(sto %>% select(strategy, shape, D, median_T, share_eliminated_50y), n = 60)
