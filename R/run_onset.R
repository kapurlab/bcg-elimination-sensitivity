# Speed of onset: the interval from birth to effective immunity, and calf
# exposure relative to adults, against time to elimination.
# Run from the project root: Rscript R/run_onset.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R")
source("R/model.R")
source("R/revaccination.R")
source("R/onset.R")
herds <- load_herds("data")
qR0 <- quantile(herds$R0, c(0.25, 0.5, 0.75, 0.9))
rep_herds <- data.frame(label = sprintf(c("R0 median (%.1f)", "R0 75th pct (%.1f)"), qR0[2:3]),
                        R0 = as.numeric(qR0[2:3]))

tau_grid <- c(0, 15, 30, 45, 60, 90, 120, 180)        # days from birth to protection
m_grid <- c(1, 3, 10)                                 # calf exposure relative to adults
regimes <- data.frame(
  regime = c("Lifelong protection, calves at birth (paper)",
             "18-month protection, annual boost of all uninfected animals"),
  D = c(Inf, 1.5), k = c(1, 20), revacc = c(FALSE, TRUE), redose_V = c(FALSE, TRUE))
design <- merge(merge(regimes, data.frame(tau_days = tau_grid)), data.frame(m = m_grid))

res <- bind_rows(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, ]
  sim <- simulate_onset(rep_herds$R0, rep(100, nrow(rep_herds)), BASE$e_s, BASE$e_i, D = d$D,
                        tau_days = d$tau_days, m = d$m, k = d$k, revacc = d$revacc,
                        redose_V = d$redose_V, horizon = 200, dt = 0.25)
  i5 <- which.min(abs(sim$times - 5))
  data.frame(d, herd = rep_herds$label, R0 = rep_herds$R0,
             T_elim = vapply(seq_len(nrow(rep_herds)), function(j)
               crossing_time(sim$times, sim$prev[, j], 1e-3), numeric(1)),
             share_calves_infected_before_onset_first5y = sim$share_infected_before_onset[i5, ])
}))
write.csv(res, "output/tables/onset.csv", row.names = FALSE)

# The whole 57-herd population, 50-year horizon, paper regime
pop <- bind_rows(lapply(seq_len(nrow(design[design$k == 1, ])), function(i) {
  d <- design[design$k == 1, ][i, ]
  sim <- simulate_onset(herds$R0, herds$herdsize, BASE$e_s, BASE$e_i, D = d$D,
                        tau_days = d$tau_days, m = d$m, k = 1, horizon = 50, dt = 0.5)
  Th <- vapply(seq_len(nrow(herds)), function(j) crossing_time(sim$times, sim$prev[, j], 1e-3), numeric(1))
  data.frame(d, share_herds_free_50y = mean(Th <= 50),
             pop_prev_50y = sum(sim$infected[nrow(sim$infected), ]) / sum(herds$herdsize))
}))
write.csv(pop, "output/tables/onset_population.csv", row.names = FALSE)

# --- figure -------------------------------------------------------------------
m_pal <- c("1" = "#7FB3D5", "3" = "#F5B7B1", "10" = "#BB8FCE")
rr <- res %>% mutate(m = factor(m), T_plot = ifelse(is.finite(T_elim), T_elim, NA),
                     regime = factor(regime, levels = regimes$regime))
g12 <- ggplot(rr, aes(tau_days, T_plot, colour = m)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_colour_manual(values = m_pal, name = "Calf exposure relative to adults") +
  scale_x_continuous(breaks = tau_grid) +
  facet_grid(herd ~ regime) +
  labs(x = "Days from birth to effective immunity", y = "Years to herd prevalence < 0.1%",
       title = "Speed of onset: the unprotected window after birth",
       subtitle = "e_s 0.58, e_i 0.74, all calves dosed. Window covers delay to vaccination plus delay to onset. Gaps: not eliminated within 200 years.") +
  theme_pub()
ggsave("output/figs/fig12_onset_delay.png", g12, width = 11, height = 7, dpi = 200, bg = "white")

rc <- res %>% filter(k == 1) %>% mutate(m = factor(m))
g13 <- ggplot(rc, aes(tau_days, share_calves_infected_before_onset_first5y, colour = m)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_colour_manual(values = m_pal, name = "Calf exposure relative to adults") +
  scale_x_continuous(breaks = tau_grid) + scale_y_continuous(labels = scales::percent) +
  facet_wrap(~herd) +
  labs(x = "Days from birth to effective immunity",
       y = "Vaccinated calves infected before onset, first 5 years of the programme",
       title = "Infections acquired in the window, paper regime") +
  theme_pub()
ggsave("output/figs/fig13_onset_window_infections.png", g13, width = 10, height = 5, dpi = 200, bg = "white")

cat("Onset done\n")
print(res %>% filter(grepl("median", herd)) %>% mutate(T_elim = round(T_elim, 1)) %>%
        select(regime, m, tau_days, T_elim) %>% pivot_wider(names_from = tau_days, values_from = T_elim), n = 20, width = 200)
