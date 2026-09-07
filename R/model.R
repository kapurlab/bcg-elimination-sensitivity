# Deterministic (mean-field) version of the within-herd transmission model of
# Fromsa et al. 2024 (Science 383:eadl3962), extended with waning of vaccine
# protection and partial calfhood coverage.
#
# Compartments follow the authors' SimInf model (transmission_model/script/
# transmission.R in the BCGCrossover repository): S susceptible, I infected,
# V vaccinated, IV infected vaccinated. Herd size is constant, births balance
# deaths, and infected animals are never removed or recovered. Two things are
# added here:
#   1. Waning: V -> S at rate w = 1/D, where D is the mean duration of
#      protection in years. D = Inf reproduces the authors' assumption of
#      lifelong protection. By default the reduced infectiousness of IV
#      animals does not wane (wane_IV = FALSE); set wane_IV = TRUE to move
#      IV -> I at the same rate.
#   2. Coverage p is the proportion of newborn calves vaccinated (the authors
#      compared p = 0 and p = 1 only).
#
# All rates are per year. Direct efficacy e_s reduces susceptibility of V,
# indirect efficacy e_i reduces infectiousness of IV, exactly as in the
# original transition list.

suppressPackageStartupMessages({
  library(deSolve)
})

U_MORT <- 0.2731918   # annual mortality rate, Demography/MortalityExp.csv

# Base case: posterior medians of the DST1 chain binomial fit, p = 1 and
# lifelong protection as assumed in the paper.
BASE <- list(e_s = 0.58, e_i = 0.74, D = Inf, p = 1)

load_herds <- function(dir = "data") {
  r0 <- read.csv(file.path(dir, "R0Estimates.csv"))
  hs <- read.csv(file.path(dir, "herdsize.csv"))
  farm <- as.integer(sub("R0.", "", names(r0), fixed = TRUE))
  data.frame(farm = farm,
             herdsize = hs$herdsize[match(farm, hs$farm)],
             R0 = apply(r0, 2, median))
}

load_efficacy_posterior <- function(dir = "data") {
  read.table(file.path(dir, "dst1_post.csv"), header = TRUE)
}

# Vaccinated fraction of the herd at the disease-free steady state.
vacc_fraction <- function(p, D, u = U_MORT) {
  w <- ifelse(is.infinite(D), 0, 1 / D)
  p * u / (u + w)
}

# Reproduction number under vaccination at the disease-free steady state.
# A vaccinated animal contributes (1 - e_s)(1 - e_i) of an unvaccinated one.
R_vacc <- function(R0, e_s, e_i, D, p, u = U_MORT) {
  f <- vacc_fraction(p, D, u)
  R0 * ((1 - f) + f * (1 - e_s) * (1 - e_i))
}

# Right-hand side for H independent herds stacked in one state vector.
rhs <- function(t, y, pr) {
  H <- pr$H
  S <- y[1:H]; I <- y[H + 1:H]; V <- y[2 * H + 1:H]; IV <- y[3 * H + 1:H]
  N <- pr$N
  lam <- pr$b0 * (I + (1 - pr$e_i) * IV) / N
  u <- pr$u; w <- pr$w; p <- pr$p; e_s <- pr$e_s
  wIV <- if (pr$wane_IV) w else 0
  dS  <- (1 - p) * u * N - lam * S - u * S + w * V
  dI  <- lam * S - u * I + wIV * IV
  dV  <- p * u * N - (1 - e_s) * lam * V - u * V - w * V
  dIV <- (1 - e_s) * lam * V - u * IV - wIV * IV
  list(c(dS, dI, dV, dIV))
}

# Simulate H herds initialised at their endemic equilibrium, then vaccinated
# from t = 0. Returns a matrix of herd prevalences (rows = times).
simulate_herds <- function(R0, N, e_s, e_i, D, p, horizon = 200, dt = 0.25,
                           u = U_MORT, wane_IV = FALSE, ...) {
  H <- length(R0)
  S0 <- N / R0; I0 <- N - S0
  I0[R0 <= 1] <- 0; S0[R0 <= 1] <- N[R0 <= 1]
  y0 <- c(S0, I0, rep(0, H), rep(0, H))
  pr <- list(H = H, N = N, b0 = R0 * u, u = u,
             w = if (is.infinite(D)) 0 else 1 / D,
             p = p, e_s = e_s, e_i = e_i, wane_IV = wane_IV)
  times <- seq(0, horizon, by = dt)
  out <- ode(y0, times, rhs, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10)
  inf <- out[, 1 + H + 1:H, drop = FALSE] + out[, 1 + 3 * H + 1:H, drop = FALSE]
  list(times = times, prev = sweep(inf, 2, N, "/"), infected = inf, N = N)
}

# First time a prevalence series drops below thr, by log-linear interpolation.
# Returns Inf if the threshold is never crossed within the horizon.
crossing_time <- function(times, prev, thr) {
  below <- which(prev < thr)
  if (length(below) == 0) return(Inf)
  k <- below[1]
  if (k == 1) return(0)
  lp <- log(prev[c(k - 1, k)])
  times[k - 1] + (log(thr) - lp[1]) / (lp[2] - lp[1]) * (times[k] - times[k - 1])
}

# Time to elimination for each herd and for the size-weighted population.
time_to_elimination <- function(R0, N, e_s, e_i, D, p, thr = 1e-3,
                                horizon = 200, ...) {
  sim <- simulate_herds(R0, N, e_s, e_i, D, p, horizon = horizon, ...)
  herd <- vapply(seq_along(R0), function(j)
    crossing_time(sim$times, sim$prev[, j], thr), numeric(1))
  pop_prev <- rowSums(sim$infected) / sum(N)
  list(herd = herd, population = crossing_time(sim$times, pop_prev, thr),
       final_prev = pop_prev[length(pop_prev)])
}

# Convenience wrapper for a single representative herd.
T_single <- function(R0, e_s, e_i, D, p, N = 100, thr = 1e-3, horizon = 200,
                     ...) {
  time_to_elimination(R0, N, e_s, e_i, D, p, thr = thr, horizon = horizon,
                      ...)$herd[1]
}

# Pastel palette used across figures.
PAL <- c(direct = "#7FB3D5", indirect = "#F5B7B1", duration = "#A9DFBF",
         coverage = "#F9E79F", R0 = "#D7BDE2")
theme_pub <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, colour = "grey70"),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
                   strip.text = ggplot2::element_text(face = "bold"),
                   legend.position = "bottom")
}
