# One-at-a-time sensitivity of time to elimination to direct efficacy,
# indirect efficacy, duration of protection and calfhood coverage.
# Run from the project root: Rscript R/run_oat.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/model.R")
herds <- load_herds("data")
post <- load_efficacy_posterior("data")

# Representative herds: quartiles and 90th percentile of the posterior-median R0
qR0 <- quantile(herds$R0, c(0.25, 0.5, 0.75, 0.9))
rep_herds <- data.frame(label = c("R0 25th pct", "R0 median", "R0 75th pct", "R0 90th pct"),
                        R0 = as.numeric(qR0))
rep_herds$label <- sprintf("%s (%.1f)", rep_herds$label, rep_herds$R0)
rep_herds$label <- factor(rep_herds$label, levels = rep_herds$label)

sweep_grid <- list(
  direct   = seq(0, 0.95, by = 0.025),
  indirect = seq(0, 0.95, by = 0.025),
  duration = c(seq(1, 10, by = 0.5), seq(11, 40, by = 1), 60, 100, Inf),
  coverage = seq(0.2, 1, by = 0.02))

set_param <- function(name, value) {
  b <- BASE
  b[[c(direct = "e_s", indirect = "e_i", duration = "D", coverage = "p")[name]]] <- value
  b
}

# --- representative herds -------------------------------------------------
oat <- bind_rows(lapply(names(sweep_grid), function(nm) {
  bind_rows(lapply(sweep_grid[[nm]], function(v) {
    b <- set_param(nm, v)
    T <- vapply(rep_herds$R0, function(r) T_single(r, b$e_s, b$e_i, b$D, b$p),
                numeric(1))
    data.frame(parameter = nm, value = v, herd = rep_herds$label,
               R0 = rep_herds$R0, T_elim = T,
               R_v = R_vacc(rep_herds$R0, b$e_s, b$e_i, b$D, b$p))
  }))
}))
write.csv(oat, "output/tables/oat_representative_herds.csv", row.names = FALSE)

# --- whole 57-herd population --------------------------------------------
oat_pop <- bind_rows(lapply(names(sweep_grid), function(nm) {
  bind_rows(lapply(sweep_grid[[nm]], function(v) {
    b <- set_param(nm, v)
    r <- time_to_elimination(herds$R0, herds$herdsize, b$e_s, b$e_i, b$D, b$p,
                             horizon = 100)
    sim <- simulate_herds(herds$R0, herds$herdsize, b$e_s, b$e_i, b$D, b$p,
                          horizon = 50, dt = 1)
    prev50 <- sum(sim$infected[nrow(sim$infected), ]) / sum(herds$herdsize)
    data.frame(parameter = nm, value = v,
               share_herds_free_50y = mean(r$herd <= 50),
               share_herds_free_100y = mean(r$herd <= 100),
               share_herds_Rv_below_1 = mean(R_vacc(herds$R0, b$e_s, b$e_i, b$D, b$p) < 1),
               pop_prev_50y = prev50,
               median_herd_T = median(r$herd))
  }))
}))
write.csv(oat_pop, "output/tables/oat_population.csv", row.names = FALSE)

# --- local elasticities at the base case ----------------------------------
# d ln T / d ln x by central differences. Duration is evaluated at D = 20 and
# 30 years because the base case (lifelong protection) has zero slope by
# construction.
elasticity <- function(R0, name, x0, h = 0.02) {
  f <- function(x) { b <- set_param(name, x); T_single(R0, b$e_s, b$e_i, b$D, b$p) }
  # coverage cannot exceed 1, so use a backward difference at p = 1
  if (name == "coverage" && x0 >= 1)
    return((log(f(x0)) - log(f(x0 * (1 - 2 * h)))) / (log(x0) - log(x0 * (1 - 2 * h))))
  (log(f(x0 * (1 + h))) - log(f(x0 * (1 - h)))) / (log(1 + h) - log(1 - h))
}
elas <- bind_rows(lapply(seq_len(nrow(rep_herds)), function(i) {
  R0 <- rep_herds$R0[i]
  data.frame(herd = rep_herds$label[i],
             parameter = c("direct", "indirect", "coverage", "duration", "duration"),
             at = c(0.58, 0.74, 1, 20, 30),
             elasticity = c(elasticity(R0, "direct", 0.58),
                            elasticity(R0, "indirect", 0.74),
                            elasticity(R0, "coverage", 1),
                            elasticity(R0, "duration", 20),
                            elasticity(R0, "duration", 30)))
}))
write.csv(elas, "output/tables/local_elasticities.csv", row.names = FALSE)

# --- where elimination is lost (R_v = 1) for each parameter -----------------
# Solve R_v(x) = 1 for the swept parameter with the others at base.
cliff <- bind_rows(lapply(rep_herds$R0, function(R0) {
  g <- function(name, x) { b <- set_param(name, x); R_vacc(R0, b$e_s, b$e_i, b$D, b$p) - 1 }
  root <- function(name, lo, hi) {
    if (g(name, lo) * g(name, hi) > 0) return(NA_real_)
    uniroot(function(x) g(name, x), c(lo, hi))$root
  }
  data.frame(R0 = R0,
             direct_min = root("direct", 0, 0.999),
             indirect_min = root("indirect", 0, 0.999),
             duration_min_years = root("duration", 0.1, 500),
             coverage_min = root("coverage", 0.01, 1))
}))
cliff$herd <- rep_herds$label
write.csv(cliff, "output/tables/elimination_thresholds.csv", row.names = FALSE)

# --- figures ---------------------------------------------------------------
lab <- c(direct = "Direct efficacy (e_s)", indirect = "Indirect efficacy (e_i)",
         duration = "Duration of protection (years)", coverage = "Calfhood coverage (p)")
base_vline <- data.frame(parameter = names(lab),
                         x = c(0.58, 0.74, NA, 1))
herd_pal <- c("#7FB3D5", "#F1948A", "#82E0AA", "#BB8FCE")

oat_plot <- oat %>% mutate(T_plot = ifelse(is.finite(T_elim), T_elim, NA),
                           parameter = factor(parameter, levels = names(lab)),
                           value_plot = ifelse(parameter == "duration" & is.infinite(value), 120, value))
p1 <- ggplot(oat_plot, aes(value_plot, T_plot, colour = herd)) +
  geom_hline(yintercept = 50, linetype = 3, colour = "grey55") +
  geom_line(linewidth = 1) +
  geom_vline(data = base_vline, aes(xintercept = x), linetype = 2, colour = "grey40") +
  scale_colour_manual(values = herd_pal, name = NULL) +
  scale_y_continuous(limits = c(0, 200)) +
  facet_wrap(~parameter, scales = "free_x", labeller = as_labeller(lab)) +
  labs(x = NULL, y = "Years to herd prevalence < 0.1%",
       title = "Time to elimination, one parameter varied at a time",
       subtitle = "Other parameters at the paper's base case (e_s 0.58, e_i 0.74, lifelong protection, 100% coverage).\nCurves stop where R_v reaches 1 and elimination is no longer possible. Duration axis: 120 = lifelong.") +
  theme_pub()
ggsave("output/figs/fig1_oat_time_to_elimination.png", p1, width = 10, height = 7.5, dpi = 200, bg = "white")

p1b <- ggplot(oat_plot, aes(value_plot, R_v, colour = herd)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey30") +
  geom_line(linewidth = 1) +
  geom_vline(data = base_vline, aes(xintercept = x), linetype = 2, colour = "grey40") +
  scale_colour_manual(values = herd_pal, name = NULL) +
  facet_wrap(~parameter, scales = "free_x", labeller = as_labeller(lab)) +
  labs(x = NULL, y = expression(R[v]~"(reproduction number under vaccination)"),
       title = "Reproduction number under vaccination, one parameter varied at a time") +
  theme_pub()
ggsave("output/figs/fig1b_oat_Rv.png", p1b, width = 10, height = 7.5, dpi = 200, bg = "white")

pop_long <- oat_pop %>%
  mutate(parameter = factor(parameter, levels = names(lab)),
         value_plot = ifelse(parameter == "duration" & is.infinite(value), 120, value)) %>%
  select(parameter, value_plot, share_herds_free_50y, share_herds_Rv_below_1, pop_prev_50y) %>%
  pivot_longer(-c(parameter, value_plot)) %>%
  mutate(name = recode(name, share_herds_free_50y = "Share of herds below 0.1% by year 50",
                       share_herds_Rv_below_1 = "Share of herds with R_v < 1",
                       pop_prev_50y = "Population animal prevalence at year 50"))
p2 <- ggplot(pop_long, aes(value_plot, value, colour = name)) +
  geom_line(linewidth = 1) +
  geom_vline(data = base_vline, aes(xintercept = x), linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c("#7FB3D5", "#F5B7B1", "#82E0AA"), name = NULL) +
  facet_wrap(~parameter, scales = "free_x", labeller = as_labeller(lab)) +
  labs(x = NULL, y = "Proportion", title = "Population outcomes across the 57 study herds") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_pub()
ggsave("output/figs/fig2_oat_population.png", p2, width = 10, height = 7.5, dpi = 200, bg = "white")

# Duration x coverage heatmap for the median herd (maroon sequential palette).
grid <- expand.grid(D = c(seq(1, 40, by = 1), Inf), p = seq(0.2, 1, by = 0.02))
grid$T <- mapply(function(D, p) T_single(qR0[2], BASE$e_s, BASE$e_i, D, p), grid$D, grid$p)
grid$T_plot <- ifelse(is.finite(grid$T), pmin(grid$T, 150), NA)
grid$D_plot <- ifelse(is.infinite(grid$D), 42, grid$D)
p3 <- ggplot(grid, aes(D_plot, p, fill = T_plot)) +
  geom_tile() +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", na.value = "grey85",
                      name = "Years to < 0.1%\n(grey: never)") +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 42), labels = c("1", "10", "20", "30", "40", "life")) +
  labs(x = "Duration of protection (years)", y = "Calfhood coverage",
       title = "Duration and coverage act through the same quantity: the vaccinated fraction of the herd",
       subtitle = sprintf("Median herd (R0 = %.2f), e_s = 0.58, e_i = 0.74", qR0[2])) +
  theme_pub()
ggsave("output/figs/fig3_duration_coverage_heatmap.png", p3, width = 9, height = 6, dpi = 200, bg = "white")

# Direct x indirect heatmap for the median herd (yellow sequential palette).
grid2 <- expand.grid(e_s = seq(0, 0.95, by = 0.025), e_i = seq(0, 0.95, by = 0.025))
grid2$T <- mapply(function(a, b) T_single(qR0[2], a, b, Inf, 1), grid2$e_s, grid2$e_i)
grid2$T_plot <- ifelse(is.finite(grid2$T), pmin(grid2$T, 150), NA)
p4 <- ggplot(grid2, aes(e_s, e_i, fill = T_plot)) +
  geom_tile() +
  annotate("point", x = 0.58, y = 0.74, shape = 21, size = 3, fill = "white") +
  scale_fill_gradient(low = "#FFF9DB", high = "#B7950B", na.value = "grey85",
                      name = "Years to < 0.1%\n(grey: never)") +
  labs(x = "Direct efficacy (e_s)", y = "Indirect efficacy (e_i)",
       title = "Direct and indirect efficacy are interchangeable through (1 - e_s)(1 - e_i)",
       subtitle = sprintf("Median herd (R0 = %.2f), lifelong protection, 100%% coverage. Point: paper's estimate.", qR0[2])) +
  theme_pub()
ggsave("output/figs/fig4_direct_indirect_heatmap.png", p4, width = 8, height = 6.5, dpi = 200, bg = "white")

cat("OAT done\n")
print(elas)
print(cliff)
