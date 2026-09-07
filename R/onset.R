# Extension of revaccination.R: a delay between vaccination at birth and the
# onset of protection. Vaccinated calves enter V0 (dosed, not yet protected)
# and move to V1 at rate 1/tau, where tau is the mean interval in years from
# birth to effective immunity. It covers both the delay to vaccination and the
# delay to onset after the dose. While in V0 a calf is fully susceptible and
# its exposure can be scaled by m relative to an adult, to represent
# exposure through milk and cohabitation. Infections acquired in V0 are
# unprotected infections (I) and are counted in C for reporting.

suppressPackageStartupMessages({library(deSolve)})

rhs_onset <- function(t, y, pr) {
  H <- pr$H; k <- pr$k; N <- pr$N
  idx <- function(j) (j - 1) * H + 1:H
  S <- y[idx(1)]; V0 <- y[idx(2)]; I <- y[idx(3)]; IV <- y[idx(4)]
  V <- matrix(y[(4 * H + 1):((4 + k) * H)], nrow = H)
  Vtot <- rowSums(V)
  lam <- pr$b0 * (I + (1 - pr$e_i) * IV) / N
  u <- pr$u; a <- pr$k_rate; g <- pr$onset_rate; m <- pr$m
  dS  <- (1 - pr$p) * u * N - lam * S - u * S + a * V[, k]
  dV0 <- pr$p * u * N * (g > 0) - g * V0 - m * lam * V0 - u * V0
  dI  <- lam * S + m * lam * V0 - u * I
  dIV <- (1 - pr$e_s) * lam * Vtot - u * IV
  dV <- matrix(0, H, k)
  for (j in 1:k) {
    inflow <- if (j == 1) (if (g > 0) g * V0 else pr$p * u * N) else a * V[, j - 1]
    dV[, j] <- inflow - a * V[, j] - (1 - pr$e_s) * lam * V[, j] - u * V[, j]
  }
  dC <- m * lam * V0                       # cumulative infections before onset
  dB <- pr$p * u * N                       # cumulative vaccinated births
  list(c(dS, dV0, dI, dIV, as.vector(dV), dC, dB))
}

pulse_onset <- function(t, y, pr) {
  H <- pr$H; k <- pr$k
  moved <- pr$cv * y[1:H]
  y[1:H] <- y[1:H] - moved
  y[4 * H + 1:H] <- y[4 * H + 1:H] + moved
  if (pr$redose_V && k > 1) for (j in 2:k) {
    idx <- (3 + j) * H + 1:H
    re <- pr$cv * y[idx]; y[idx] <- y[idx] - re; y[4 * H + 1:H] <- y[4 * H + 1:H] + re
  }
  y
}

simulate_onset <- function(R0, N, e_s, e_i, D, tau_days, m = 1, p = 1, k = 1,
                           interval = 1, cv = 1, horizon = 200, dt = 0.25,
                           u = U_MORT, revacc = FALSE, redose_V = FALSE) {
  H <- length(R0)
  S0 <- N / R0; I0 <- N - S0
  y0 <- c(S0, rep(0, H), I0, rep(0, H), rep(0, H * k), rep(0, H), rep(0, H))
  pr <- list(H = H, N = N, b0 = R0 * u, u = u, k = k,
             k_rate = if (is.infinite(D)) 0 else k / D,
             onset_rate = if (tau_days > 0) 365 / tau_days else 0, m = m,
             p = p, e_s = e_s, e_i = e_i, cv = cv, redose_V = redose_V)
  times <- seq(0, horizon, by = dt)
  if (revacc) {
    ev <- seq(interval, horizon, by = interval); times <- sort(unique(c(times, ev)))
    out <- ode(y0, times, rhs_onset, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10,
               events = list(func = pulse_onset, time = ev))
  } else out <- ode(y0, times, rhs_onset, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10)
  tt <- out[, 1]
  inf <- out[, 1 + 2 * H + 1:H, drop = FALSE] + out[, 1 + 3 * H + 1:H, drop = FALSE]
  C <- out[, 1 + (4 + k) * H + 1:H, drop = FALSE]; B <- out[, 1 + (5 + k) * H + 1:H, drop = FALSE]
  list(times = tt, prev = sweep(inf, 2, N, "/"), infected = inf, N = N,
       share_infected_before_onset = C / pmax(B, 1e-9))
}
