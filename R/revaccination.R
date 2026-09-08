# Extension of model.R: annual revaccination of every susceptible animal in
# addition to vaccination of calves at birth, with the duration of protection
# modelled as an Erlang distribution with k stages (V1 ... Vk). k = 1 is the
# exponential waning of model.R. Large k (20 here) makes the duration nearly
# fixed, so protection lasts close to exactly D years after each dose.
#
# Revaccination is a pulse: at each campaign a fraction `cv` of the S
# compartment moves to V1. Infected animals gain nothing from vaccination in
# this model, so they are not moved.

suppressPackageStartupMessages({library(deSolve)})

rhs_stage <- function(t, y, pr) {
  H <- pr$H; k <- pr$k; N <- pr$N
  idx <- function(j) (j - 1) * H + 1:H
  S <- y[idx(1)]; I <- y[idx(2)]; IV <- y[idx(3)]
  V <- matrix(y[(3 * H + 1):((3 + k) * H)], nrow = H)   # H x k
  Vtot <- rowSums(V)
  lam <- pr$b0 * (I + (1 - pr$e_i) * IV) / N
  u <- pr$u; a <- pr$k_rate                              # a = k / D
  dS  <- (1 - pr$p) * u * N - lam * S - u * S + a * V[, k]
  dI  <- lam * S - u * I
  dIV <- (1 - pr$e_s) * lam * Vtot - u * IV
  dV <- matrix(0, H, k)
  for (j in 1:k) {
    dV[, j] <- (if (j == 1) pr$p * u * N else a * V[, j - 1]) -
      a * V[, j] - (1 - pr$e_s) * lam * V[, j] - u * V[, j]
  }
  list(c(dS, dI, dIV, as.vector(dV)))
}

# Pulse: move a fraction cv of S into V1 at each campaign. With redose_V,
# animals still under protection are dosed too, which restarts their
# protection clock (all V stages collapse into V1).
pulse_event <- function(t, y, pr) {
  H <- pr$H; k <- pr$k
  moved <- pr$cv * y[1:H]
  y[1:H] <- y[1:H] - moved
  y[3 * H + 1:H] <- y[3 * H + 1:H] + moved
  if (pr$redose_V && k > 1) {
    for (j in 2:k) {
      idx <- (2 + j) * H + 1:H
      re <- pr$cv * y[idx]
      y[idx] <- y[idx] - re
      y[3 * H + 1:H] <- y[3 * H + 1:H] + re
    }
  }
  y
}

simulate_revacc <- function(R0, N, e_s, e_i, D, p = 1, k = 1, interval = 1,
                            cv = 1, horizon = 200, dt = 0.25, u = U_MORT,
                            revacc = TRUE, redose_V = FALSE) {
  H <- length(R0)
  S0 <- N / R0; I0 <- N - S0
  I0[R0 <= 1] <- 0; S0[R0 <= 1] <- N[R0 <= 1]
  y0 <- c(S0, I0, rep(0, H), rep(0, H * k))
  pr <- list(H = H, N = N, b0 = R0 * u, u = u, k = k,
             k_rate = if (is.infinite(D)) 0 else k / D,
             p = p, e_s = e_s, e_i = e_i, cv = cv, redose_V = redose_V)
  times <- seq(0, horizon, by = dt)
  if (revacc) {
    ev_times <- seq(interval, horizon, by = interval)
    times <- sort(unique(c(times, ev_times)))
    out <- ode(y0, times, rhs_stage, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10,
               events = list(func = pulse_event, time = ev_times))
  } else {
    out <- ode(y0, times, rhs_stage, pr, method = "lsoda", rtol = 1e-8, atol = 1e-10)
  }
  tt <- out[, 1]
  inf <- out[, 1 + H + 1:H, drop = FALSE] + out[, 1 + 2 * H + 1:H, drop = FALSE]
  Vtot <- out[, 1 + 3 * H + 1:(H * k), drop = FALSE]
  Vtot <- sapply(1:H, function(i) rowSums(Vtot[, i + (0:(k - 1)) * H, drop = FALSE]))
  list(times = tt, prev = sweep(inf, 2, N, "/"), infected = inf, N = N,
       vfrac = sweep(Vtot, 2, N, "/"))
}

# Time-averaged vaccinated fraction of a disease-free herd, over the last
# campaign cycle of a 30-year run, and the corresponding average R_v.
vacc_fraction_revacc <- function(D, k = 1, interval = 1, cv = 1, p = 1, u = U_MORT,
                                 revacc = TRUE, redose_V = FALSE) {
  sim <- simulate_revacc(R0 = 1, N = 1000, e_s = 0, e_i = 0, D = D, p = p, k = k,
                         interval = interval, cv = cv, horizon = 30, dt = 0.05,
                         revacc = revacc, redose_V = redose_V)
  # R0 = 1 with a disease-free start keeps the herd uninfected
  sel <- sim$times > 30 - interval
  mean(sim$vfrac[sel, 1])
}

# ---- stochastic version (SimInf) ------------------------------------------
# Same structure as run_stochastic.R with k protection stages and annual
# campaigns implemented as internal transfer events. Returns extinction
# times in years, one per replicate herd.
stoch_revacc <- function(R0, N, e_s, e_i, D, k = 1, p = 1, interval = 1, cv = 1,
                         reps = 300, years = 200, u = U_MORT, revacc = TRUE,
                         redose_V = FALSE) {
  suppressPackageStartupMessages(library(SimInf))
  Vn <- paste0("V", 1:k)
  comp <- c("S", Vn, "I", "IV")
  Vsum <- paste(Vn, collapse = "+")
  Ntot <- paste(comp, collapse = "+")
  inf_term <- sprintf("b0*(I + (1-e_i)*IV)/(%s)", Ntot)
  tr <- c(sprintf("@ -> (1-p)*br*N0 -> S"),
          sprintf("@ -> p*br*N0 -> V1"),
          sprintf("S -> S > 0 ? %s*S : 0 -> I", inf_term),
          sprintf("S -> S > 0 ? S*mu : 0 -> @"),
          sprintf("I -> I > 0 ? I*mu : 0 -> @"),
          sprintf("IV -> IV > 0 ? IV*mu : 0 -> @"))
  for (j in 1:k) {
    nxt <- if (j < k) Vn[j + 1] else "S"
    tr <- c(tr,
            sprintf("%s -> %s > 0 ? (1-e_s)*%s*%s : 0 -> IV", Vn[j], Vn[j], inf_term, Vn[j]),
            sprintf("%s -> %s > 0 ? %s*mu : 0 -> @", Vn[j], Vn[j], Vn[j]),
            sprintf("%s -> %s > 0 ? %s*w : 0 -> %s", Vn[j], Vn[j], Vn[j], nxt))
  }
  mu <- u / 365
  R0 <- rep(R0, reps); N <- rep(N, reps)
  S0 <- round(N / R0); I0 <- N - S0
  u0 <- as.data.frame(matrix(0, length(R0), length(comp), dimnames = list(NULL, comp)))
  u0$S <- S0; u0$I <- I0
  ldata <- data.frame(N0 = N, b0 = R0 * mu, br = mu, mu = mu)
  gdata <- c(p = p, e_s = e_s, e_i = e_i, w = if (is.infinite(D)) 0 else k / (365 * D))
  # E: column 1 selects S (and the V stages when redosing); N: shift to V1
  E <- matrix(0, length(comp), 1, dimnames = list(comp, "1"))
  E["S", 1] <- 1
  if (redose_V && k > 1) E[Vn[-1], 1] <- 1
  Nm <- matrix(0, length(comp), 1, dimnames = list(comp, "1"))
  Nm["S", 1] <- 1                                   # S -> V1 is one step forward
  if (redose_V && k > 1) Nm[Vn[-1], 1] <- -(1:(k - 1))  # Vj -> V1
  events <- NULL
  if (revacc) {
    ev_days <- as.integer(round(seq(365 * interval, 365 * years, by = 365 * interval)))
    events <- data.frame(event = "intTrans",
                         time = rep(ev_days, each = length(R0)),
                         node = as.integer(rep(seq_along(R0), times = length(ev_days))),
                         dest = 0L, n = 0L, proportion = cv, select = 1L, shift = 1L)
  }
  tspan <- seq(1, 365 * years, by = 30)
  m <- mparse(transitions = tr, compartments = comp, ldata = ldata, gdata = gdata,
              u0 = u0, tspan = tspan, events = events, E = E, N = Nm)
  x <- trajectory(run(m))
  x$inf <- x$I + x$IV
  x$rep <- x$node
  x %>% dplyr::group_by(rep) %>%
    dplyr::summarise(T = if (any(inf == 0)) min(time[inf == 0]) / 365 else Inf) %>%
    dplyr::pull(T)
}
