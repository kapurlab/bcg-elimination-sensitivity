# Extension of boost.R: protection that varies with time since dose.
#
# Time since dose is tracked by a chain of k stages of equal mean length
# (stage_len years each, so the chain spans k * stage_len years). Stage j
# carries a relative protection prof[j] in [0, 1]; the susceptibility
# reduction in that stage is e_s * prof[j]. An animal leaving the last stage
# becomes an unprotected susceptible. A birth dose and a booster dose can have
# different profiles. Campaign pulses move every uninfected animal to stage 1
# of the booster chain. An animal infected in a stage with relative
# protection prof enters IV (indirect effect e_i) with probability prof and
# I otherwise; this is an assumption about partial protection.
#
# Two extra ingredients used by the later scripts:
#   nonresp : fraction of calves that never respond to any dose (persistent
#             non-responders). They are born into U, are never dosed, and are
#             otherwise identical to S.
#   imports : infected animals entering the herd per year (constant N: the
#             same number of births is withheld).

suppressPackageStartupMessages({library(deSolve)})

# Rise, plateau, decline profile evaluated at the stage midpoints.
#   onset_days : days from dose to full protection (linear rise)
#   plateau_yr : years from dose at which decline begins
#   zero_yr    : years from dose at which protection reaches zero (linear)
make_profile <- function(k, stage_len, onset_days = 30, plateau_yr = 1, zero_yr = Inf) {
  t <- (seq_len(k) - 0.5) * stage_len
  rise <- pmin(1, t / (onset_days / 365))
  fall <- if (is.infinite(zero_yr)) rep(1, k) else
    pmax(0, pmin(1, (zero_yr - t) / (zero_yr - plateau_yr)))
  pmin(rise, fall)
}

# Profile from a trial hazard curve: HR(t) is the hazard ratio of vaccinated
# to unvaccinated at time t since dose; relative protection is
# (1 - HR) / e_s_max so that the plateau maps to prof = 1.
profile_from_hazard <- function(k, stage_len, t_hr, HR, e_s_max) {
  t <- (seq_len(k) - 0.5) * stage_len
  hr <- approx(t_hr, HR, xout = t, rule = 2)$y
  pmax(0, pmin(1, (1 - hr) / e_s_max))
}

rhs_profile <- function(t, y, pr) {
  H <- pr$H; k <- pr$k; N <- pr$N
  S <- y[1:H]; U <- y[H + 1:H]; I <- y[2 * H + 1:H]; IV <- y[3 * H + 1:H]
  P <- matrix(y[4 * H + 1:(H * k)], nrow = H)
  B <- matrix(y[(4 + k) * H + 1:(H * k)], nrow = H)
  lam <- pr$b0 * (I + (1 - pr$e_i) * IV) / N
  u <- pr$u; a <- pr$a; a_out <- if (pr$absorb) 0 else a   # last stage absorbing for lifelong profiles
  infP <- sweep(P, 2, 1 - pr$e_s * pr$prof_p, "*") * lam   # infections per stage
  infB <- sweep(B, 2, 1 - pr$e_s * pr$prof_b, "*") * lam
  births <- u * N - pr$imp
  dS  <- (1 - pr$p) * (1 - pr$q) * births - lam * S - u * S + a_out * P[, k] + a_out * B[, k]
  dU  <- pr$q * births - lam * U - u * U
  # an animal infected in a stage with relative protection prof enters IV
  # (reduced infectiousness) with probability prof and I otherwise, so an
  # animal whose protection has fully lapsed is an ordinary infection
  dI  <- lam * S + lam * U - u * I + pr$imp +
    rowSums(sweep(infP, 2, 1 - pr$prof_p, "*")) + rowSums(sweep(infB, 2, 1 - pr$prof_b, "*"))
  dIV <- rowSums(sweep(infP, 2, pr$prof_p, "*")) + rowSums(sweep(infB, 2, pr$prof_b, "*")) - u * IV
  dP <- matrix(0, H, k); dB <- matrix(0, H, k)
  for (j in 1:k) {
    aj <- if (j == k) a_out else a
    dP[, j] <- (if (j == 1) pr$p * (1 - pr$q) * births else a * P[, j - 1]) - aj * P[, j] - infP[, j] - u * P[, j]
    dB[, j] <- (if (j == 1) 0 else a * B[, j - 1]) - aj * B[, j] - infB[, j] - u * B[, j]
  }
  list(c(dS, dU, dI, dIV, as.vector(dP), as.vector(dB)))
}

pulse_profile <- function(t, y, pr) {
  H <- pr$H; k <- pr$k; cv <- pr$cv
  b1 <- (4 + k) * H + 1:H
  take <- function(idx) { m <- cv * y[idx]; y[idx] <<- y[idx] - m; y[b1] <<- y[b1] + m }
  take(1:H)
  for (j in 1:k) take(4 * H + (j - 1) * H + 1:H)
  if (k > 1) for (j in 2:k) take((4 + k) * H + (j - 1) * H + 1:H)
  y
}

simulate_profile <- function(R0, N, e_s, e_i, prof_p, prof_b = prof_p, stage_len = 0.1,
                             interval = NA, cv = 1, p = 1, q = 0, imports = 0,
                             absorb = FALSE, horizon = 200, dt = 0.25, u = U_MORT) {
  H <- length(R0); k <- length(prof_p)
  stopifnot(length(prof_b) == k)
  S0 <- N / R0; I0 <- N - S0
  y0 <- c(S0, rep(0, H), I0, rep(0, H), rep(0, 2 * H * k))
  pr <- list(H = H, N = N, b0 = R0 * u, u = u, k = k, a = 1 / stage_len,
             p = p, q = q, e_s = e_s, e_i = e_i, cv = cv, imp = imports,
             prof_p = prof_p, prof_b = prof_b, absorb = absorb)
  times <- seq(0, horizon, by = dt)
  if (!is.na(interval)) {
    ev <- seq(interval, horizon, by = interval); times <- sort(unique(c(times, ev)))
    out <- ode(y0, times, rhs_profile, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10,
               events = list(func = pulse_profile, time = ev))
  } else out <- ode(y0, times, rhs_profile, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10)
  inf <- out[, 1 + 2 * H + 1:H, drop = FALSE] + out[, 1 + 3 * H + 1:H, drop = FALSE]
  list(times = out[, 1], prev = sweep(inf, 2, N, "/"), infected = inf, N = N)
}
