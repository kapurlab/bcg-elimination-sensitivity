# Stochastic check using the authors' SimInf transition list with a waning
# event added. Herds are independent (no movement network). Elimination is
# the first day on which no infected animal remains.
# Run from the project root: Rscript R/run_stochastic.R

suppressPackageStartupMessages({library(SimInf); library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R")
source("R/model.R")
set.seed(20260907)
herds <- load_herds("data")
R0_med <- median(herds$R0)

compartments <- c("S", "I", "V", "IV")
# Original eight transitions from transmission_model/script/transmission.R,
# with p, e_s, e_i and w as global data. Rates are per day, as in the original.
transitions <- c(
  "@ -> (1-p)*br*N0 -> S",
  "S -> S > 0 ? (b0*(I + (1-e_i)*IV)*S)/(S+I+V+IV) : 0 -> I",
  "@ -> p*br*N0 -> V",
  "V -> V > 0 ? ((1-e_s)*b0*(I + (1-e_i)*IV)*V)/(S+I+V+IV) : 0 -> IV",
  "S -> S > 0 ? S*k : 0 -> @",
  "I -> I > 0 ? I*k : 0 -> @",
  "V -> V > 0 ? V*k : 0 -> @",
  "IV -> IV > 0 ? IV*k : 0 -> @",
  "V -> V > 0 ? V*w : 0 -> S")            # added: waning of protection

run_stoch <- function(R0, N, e_s, e_i, D, p, reps, years = 200, u = U_MORT) {
  # Each replicate is a copy of every herd; SimInf simulates all nodes in C.
  R0 <- rep(R0, reps); N <- rep(N, reps)
  k <- u / 365
  S0 <- round(N / R0); I0 <- N - S0
  u0 <- data.frame(S = S0, I = I0, V = 0, IV = 0)
  ldata <- data.frame(N0 = N, b0 = R0 * k, br = k, k = k)
  gdata <- c(p = p, e_s = e_s, e_i = e_i, w = if (is.infinite(D)) 0 else 1 / (365 * D))
  tspan <- seq(1, 365 * years, by = 30)
  m <- mparse(transitions = transitions, compartments = compartments,
              ldata = ldata, gdata = gdata, u0 = u0, tspan = tspan)
  tr <- trajectory(run(m))
  tr$inf <- tr$I + tr$IV
  tr$rep <- (tr$node - 1) %/% (length(R0) / reps) + 1
  tr
}

# Extinction time per replicate for a set of herds, in years.
extinction_time <- function(tr, years) {
  tr %>% group_by(rep, time) %>% summarise(inf = sum(inf), .groups = "drop") %>%
    group_by(rep) %>% summarise(T = if (any(inf == 0)) min(time[inf == 0]) / 365 else Inf)
}

# --- OAT sweeps for the median herd (N = 44), 300 replicates each ---------
grid <- list(direct = seq(0, 0.9, by = 0.1), indirect = seq(0, 0.9, by = 0.1),
             duration = c(2, 3, 5, 7, 10, 15, 20, 30, Inf), coverage = seq(0.4, 1, by = 0.1))
set_param <- function(name, value) {
  b <- BASE
  b[[c(direct = "e_s", indirect = "e_i", duration = "D", coverage = "p")[name]]] <- value
  b
}
res <- bind_rows(lapply(names(grid), function(nm) bind_rows(lapply(grid[[nm]], function(v) {
  b <- set_param(nm, v)
  tr <- run_stoch(R0_med, 44, b$e_s, b$e_i, b$D, b$p, reps = 300)
  et <- extinction_time(tr)
  det <- T_single(R0_med, b$e_s, b$e_i, b$D, b$p, N = 44)
  data.frame(parameter = nm, value = v,
             median_T = median(et$T), q25 = quantile(et$T, 0.25), q75 = quantile(et$T, 0.75),
             share_eliminated_50y = mean(et$T <= 50), share_eliminated_200y = mean(is.finite(et$T)),
             deterministic_T = det)
}))))
write.csv(res, "output/tables/stochastic_oat_median_herd.csv", row.names = FALSE)

# --- full 57-herd population at the base case, 100 replicates --------------
tr <- run_stoch(herds$R0, herds$herdsize, BASE$e_s, BASE$e_i, BASE$D, BASE$p, reps = 100)
pop <- extinction_time(tr)
herd_ext <- tr %>% mutate(herd = (node - 1) %% 57 + 1) %>% group_by(rep, herd) %>%
  summarise(T = if (any(inf == 0)) min(time[inf == 0]) / 365 else Inf, .groups = "drop") %>%
  group_by(herd) %>% summarise(median_T = median(T), share_50y = mean(T <= 50), share_200y = mean(is.finite(T)))
herd_ext$R0 <- herds$R0; herd_ext$herdsize <- herds$herdsize
write.csv(herd_ext, "output/tables/stochastic_population_base_by_herd.csv", row.names = FALSE)
pop_summary <- data.frame(share_reps_population_free_50y = mean(pop$T <= 50),
                          share_reps_population_free_100y = mean(pop$T <= 100),
                          share_reps_population_free_200y = mean(is.finite(pop$T)),
                          median_population_T = median(pop$T))
write.csv(pop_summary, "output/tables/stochastic_population_base_summary.csv", row.names = FALSE)

# --- figure ---------------------------------------------------------------------
lab <- c(direct = "Direct efficacy (e_s)", indirect = "Indirect efficacy (e_i)",
         duration = "Duration of protection (years)", coverage = "Calfhood coverage (p)")
rp <- res %>% mutate(parameter = factor(parameter, levels = names(lab)),
                     value_plot = ifelse(parameter == "duration" & is.infinite(value), 40, value),
                     median_T = ifelse(is.finite(median_T), median_T, NA),
                     q75 = ifelse(is.finite(q75), q75, NA),
                     deterministic_T = ifelse(is.finite(deterministic_T), deterministic_T, NA))
g <- ggplot(rp, aes(value_plot)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#7FB3D5", alpha = 0.35) +
  geom_line(aes(y = median_T, colour = "Stochastic median (IQR shaded)"), linewidth = 1) +
  geom_line(aes(y = deterministic_T, colour = "Deterministic, prevalence < 0.1%"), linewidth = 1, linetype = 2) +
  scale_colour_manual(values = c("#2E86C1", "#C0392B"), name = NULL) +
  facet_wrap(~parameter, scales = "free_x", labeller = as_labeller(lab)) +
  labs(x = NULL, y = "Years to zero infected animals",
       title = sprintf("Stochastic check: SimInf model with waning added, median herd (R0 = %.2f, N = 44), 300 replicates", R0_med),
       subtitle = "Duration axis: 40 = lifelong. Missing points: elimination not reached within 200 years in half of replicates.") +
  theme_pub()
ggsave("output/figs/fig8_stochastic_check.png", g, width = 10, height = 7.5, dpi = 200, bg = "white")
cat("Stochastic done\n"); print(res); print(pop_summary)
