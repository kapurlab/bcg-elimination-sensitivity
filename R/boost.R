# Extension of revaccination.R: a birth dose with its own duration D1 and
# recurring whole-herd booster campaigns with their own duration D2.
# Protection after each kind of dose is Erlang with k stages (near-fixed for
# k = 20). At every campaign a fraction cv of the uninfected herd, protected
# or not, receives a booster and enters the first booster stage. Dosing an
# animal whose protection is still active is treated as neutral. Infected
# animals gain nothing from a dose.
#
# State per herd: S, I, IV, P1..Pk (birth-dose protection), B1..Bk (booster
# protection).

suppressPackageStartupMessages({library(deSolve)})

rhs_boost <- function(t, y, pr) {
  H <- pr$H; k <- pr$k; N <- pr$N
  S <- y[1:H]; I <- y[H + 1:H]; IV <- y[2 * H + 1:H]
  P <- matrix(y[3 * H + 1:(H * k)], nrow = H)
  B <- matrix(y[(3 + k) * H + 1:(H * k)], nrow = H)
  lam <- pr$b0 * (I + (1 - pr$e_i) * IV) / N
  u <- pr$u; a1 <- pr$a1; a2 <- pr$a2; es <- 1 - pr$e_s
  Vtot <- rowSums(P) + rowSums(B)
  dS  <- (1 - pr$p) * u * N - lam * S - u * S + a1 * P[, k] + a2 * B[, k]
  dI  <- lam * S - u * I
  dIV <- es * lam * Vtot - u * IV
  dP <- matrix(0, H, k); dB <- matrix(0, H, k)
  for (j in 1:k) {
    dP[, j] <- (if (j == 1) pr$p * u * N else a1 * P[, j - 1]) - a1 * P[, j] - es * lam * P[, j] - u * P[, j]
    dB[, j] <- (if (j == 1) 0 else a2 * B[, j - 1]) - a2 * B[, j] - es * lam * B[, j] - u * B[, j]
  }
  list(c(dS, dI, dIV, as.vector(dP), as.vector(dB)))
}

pulse_boost <- function(t, y, pr) {
  H <- pr$H; k <- pr$k; cv <- pr$cv
  b1 <- (3 + k) * H + 1:H
  take <- function(idx) { m <- cv * y[idx]; y[idx] <<- y[idx] - m; y[b1] <<- y[b1] + m }
  take(1:H)                                        # S
  for (j in 1:k) take(3 * H + (j - 1) * H + 1:H)   # birth-dose stages
  if (k > 1) for (j in 2:k) take((3 + k) * H + (j - 1) * H + 1:H)   # booster stages 2..k
  y
}

simulate_boost <- function(R0, N, e_s, e_i, D1, D2, interval = NA, cv = 1, p = 1,
                           k = 20, horizon = 200, dt = 0.25, u = U_MORT) {
  H <- length(R0)
  S0 <- N / R0; I0 <- N - S0
  y0 <- c(S0, I0, rep(0, H), rep(0, 2 * H * k))
  pr <- list(H = H, N = N, b0 = R0 * u, u = u, k = k, p = p, e_s = e_s, e_i = e_i, cv = cv,
             a1 = if (is.infinite(D1)) 0 else k / D1, a2 = if (is.infinite(D2)) 0 else k / D2)
  times <- seq(0, horizon, by = dt)
  if (!is.na(interval)) {
    ev <- seq(interval, horizon, by = interval); times <- sort(unique(c(times, ev)))
    out <- ode(y0, times, rhs_boost, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10,
               events = list(func = pulse_boost, time = ev))
  } else out <- ode(y0, times, rhs_boost, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10)
  inf <- out[, 1 + H + 1:H, drop = FALSE] + out[, 1 + 2 * H + 1:H, drop = FALSE]
  list(times = out[, 1], prev = sweep(inf, 2, N, "/"), infected = inf, N = N)
}
