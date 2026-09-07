# Global sensitivity analysis: Latin hypercube sampling with partial rank
# correlation coefficients, and variance-based Sobol indices.
# Run from the project root: Rscript R/run_global.R

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(lhs); library(sensitivity)
})
source("R/root.R")
source("R/model.R")
set.seed(20260907)
herds <- load_herds("data")
R0_med <- median(herds$R0)
HORIZON <- 200

# Input ranges. Efficacies span the plausible range for BCG in cattle.
# Duration is log-uniform because a 2-year and a 40-year vaccine differ by
# orders of magnitude in waning rate. Coverage covers partial programmes.
ranges <- list(e_s = c(0.2, 0.9), e_i = c(0, 0.95), D = c(2, 40), p = c(0.3, 1))
to_params <- function(X) {
  data.frame(e_s = ranges$e_s[1] + X[, 1] * diff(ranges$e_s),
             e_i = ranges$e_i[1] + X[, 2] * diff(ranges$e_i),
             D = exp(log(ranges$D[1]) + X[, 3] * diff(log(ranges$D))),
             p = ranges$p[1] + X[, 4] * diff(ranges$p))
}

# Three outcomes for the median herd: time to elimination (censored at the
# horizon), R_v, and log10 prevalence at year 50.
outcomes <- function(P, R0 = R0_med, N = 44) {
  t(mapply(function(e_s, e_i, D, p) {
    sim <- simulate_herds(R0, N, e_s, e_i, D, p, horizon = HORIZON, dt = 0.5)
    prev <- sim$prev[, 1]
    Tel <- crossing_time(sim$times, prev, 1e-3)
    c(T_elim = min(Tel, HORIZON),
      R_v = R_vacc(R0, e_s, e_i, D, p),
      log10_prev50 = log10(prev[which.min(abs(sim$times - 50))]))
  }, P$e_s, P$e_i, P$D, P$p))
}

# --- LHS + PRCC --------------------------------------------------------------
n_lhs <- 3000
P <- to_params(randomLHS(n_lhs, 4))
Y <- outcomes(P)
lhs_out <- cbind(P, Y)
write.csv(lhs_out, "output/tables/lhs_samples.csv", row.names = FALSE)

prcc_tab <- bind_rows(lapply(colnames(Y), function(y) {
  pc <- pcc(P, Y[, y], rank = TRUE, nboot = 300)
  data.frame(outcome = y, parameter = rownames(pc$PRCC),
             prcc = pc$PRCC[, 1], lo = pc$PRCC[, 4], hi = pc$PRCC[, 5])
}))
write.csv(prcc_tab, "output/tables/prcc.csv", row.names = FALSE)

# --- Sobol (Jansen estimator) for the median herd -------------------------
n_sob <- 8000
X1 <- data.frame(randomLHS(n_sob, 4)); X2 <- data.frame(randomLHS(n_sob, 4))
sob_model <- function(X) outcomes(to_params(as.matrix(X)))[, "T_elim"]
sob <- soboljansen(model = sob_model, X1 = X1, X2 = X2, nboot = 100)
sob_tab <- data.frame(parameter = c("e_s", "e_i", "D", "p"),
                      first = sob$S[, 1], first_lo = sob$S[, 4], first_hi = sob$S[, 5],
                      total = sob$T[, 1], total_lo = sob$T[, 4], total_hi = sob$T[, 5],
                      outcome = "T_elim (median herd)")

# Sobol on log10 prevalence at 50 years (uncensored outcome)
sob_model2 <- function(X) outcomes(to_params(as.matrix(X)))[, "log10_prev50"]
sob2 <- soboljansen(model = sob_model2, X1 = X1, X2 = X2, nboot = 100)
sob_tab2 <- data.frame(parameter = c("e_s", "e_i", "D", "p"),
                       first = sob2$S[, 1], first_lo = sob2$S[, 4], first_hi = sob2$S[, 5],
                       total = sob2$T[, 1], total_lo = sob2$T[, 4], total_hi = sob2$T[, 5],
                       outcome = "log10 prevalence at year 50 (median herd)")

# Sobol with herd R0 as a fifth factor, log-uniform over the observed range
# of posterior-median herd R0.
R0_range <- range(herds$R0)
outcomes5 <- function(X) {
  P <- to_params(X[, 1:4])
  R0 <- exp(log(R0_range[1]) + X[, 5] * diff(log(R0_range)))
  mapply(function(e_s, e_i, D, p, r) {
    sim <- simulate_herds(r, 44, e_s, e_i, D, p, horizon = HORIZON, dt = 0.5)
    min(crossing_time(sim$times, sim$prev[, 1], 1e-3), HORIZON)
  }, P$e_s, P$e_i, P$D, P$p, R0)
}
X1b <- data.frame(randomLHS(n_sob, 5)); X2b <- data.frame(randomLHS(n_sob, 5))
sob3 <- soboljansen(model = function(X) outcomes5(as.matrix(X)), X1 = X1b, X2 = X2b, nboot = 100)
sob_tab3 <- data.frame(parameter = c("e_s", "e_i", "D", "p", "R0"),
                       first = sob3$S[, 1], first_lo = sob3$S[, 4], first_hi = sob3$S[, 5],
                       total = sob3$T[, 1], total_lo = sob3$T[, 4], total_hi = sob3$T[, 5],
                       outcome = "T_elim with herd R0 as a factor")
sobol_all <- bind_rows(sob_tab, sob_tab2, sob_tab3)
write.csv(sobol_all, "output/tables/sobol.csv", row.names = FALSE)

# --- figures -----------------------------------------------------------------
plab <- c(e_s = "Direct efficacy", e_i = "Indirect efficacy", D = "Duration", p = "Coverage", R0 = "Herd R0")
pal <- c("Direct efficacy" = "#7FB3D5", "Indirect efficacy" = "#F5B7B1", "Duration" = "#A9DFBF",
         "Coverage" = "#F9E79F", "Herd R0" = "#D7BDE2")
olab <- c(T_elim = "Time to elimination (censored at 200 y)", R_v = "R_v", log10_prev50 = "log10 prevalence at year 50")

pp <- prcc_tab %>% mutate(parameter = plab[parameter], outcome = olab[outcome])
g1 <- ggplot(pp, aes(x = reorder(parameter, abs(prcc)), y = prcc, fill = parameter)) +
  geom_col(colour = "grey40", width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_hline(yintercept = 0) + coord_flip() +
  scale_fill_manual(values = pal, guide = "none") +
  facet_wrap(~outcome) +
  labs(x = NULL, y = "Partial rank correlation coefficient (95% bootstrap CI)",
       title = "Global sensitivity: PRCC from 3000 Latin hypercube samples",
       subtitle = sprintf("Median herd (R0 = %.2f). Ranges: e_s 0.2-0.9, e_i 0-0.95, duration 2-40 y (log-uniform), coverage 0.3-1.", R0_med)) +
  theme_pub()
ggsave("output/figs/fig5_prcc.png", g1, width = 11, height = 4.8, dpi = 200, bg = "white")

ss <- sobol_all %>% mutate(parameter = plab[parameter]) %>%
  pivot_longer(c(first, total), names_to = "index", values_to = "value") %>%
  mutate(lo = ifelse(index == "first", first_lo, total_lo),
         hi = ifelse(index == "first", first_hi, total_hi),
         index = recode(index, first = "First order", total = "Total effect"))
g2 <- ggplot(ss, aes(x = parameter, y = value, fill = parameter, alpha = index)) +
  geom_col(position = position_dodge(width = 0.8), colour = "grey40", width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_alpha_manual(values = c(0.55, 1), name = NULL) +
  facet_wrap(~outcome, scales = "free_x") +
  labs(x = NULL, y = "Sobol index", title = "Variance decomposition (Sobol, Jansen estimator, n = 8000 x (k + 2))") +
  theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("output/figs/fig6_sobol.png", g2, width = 12, height = 5, dpi = 200, bg = "white")

# Scatter of LHS outcomes against each input, to show shape not just rank.
sc <- lhs_out %>% mutate(T_plot = ifelse(T_elim >= HORIZON, NA, T_elim)) %>%
  pivot_longer(c(e_s, e_i, D, p), names_to = "parameter") %>%
  mutate(parameter = factor(plab[parameter], levels = plab[1:4]))
g3 <- ggplot(sc, aes(value, T_plot, colour = R_v < 1)) +
  geom_point(alpha = 0.35, size = 0.9) +
  scale_colour_manual(values = c(`TRUE` = "#7FB3D5", `FALSE` = "#F1948A"),
                      labels = c(`TRUE` = "R_v < 1", `FALSE` = "R_v >= 1 (never eliminates)"), name = NULL) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "Years to < 0.1% (blank: not within 200 y)",
       title = "Latin hypercube samples: time to elimination against each input, median herd") +
  theme_pub()
ggsave("output/figs/fig7_lhs_scatter.png", g3, width = 10, height = 7, dpi = 200, bg = "white")

cat("Global done\n"); print(prcc_tab); print(sobol_all)
