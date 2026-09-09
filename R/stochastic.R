# General stochastic (SimInf) single-herd scenario runner used for the memo
# figures. Same transition structure as the authors' model with:
#   k-stage near-fixed protection (V1..Vk), waning to S
#   coverage p of calves dosed at birth
#   an optional pre-immunity window: dosed calves enter V0 and become
#     protected at rate 365/tau_days; while in V0 they are fully susceptible
#     with exposure m_calf times an adult's
#   whole-herd booster campaigns every `interval` years that move S and
#     V2..Vk to V1 (V0 calves are not re-dosed)
# Returns the first year each replicate herd has no infected animal (Inf if
# never). With detail = TRUE it returns a data frame that also carries the
# mean within-herd prevalence over the first 20 and 50 years, a continuous
# outcome that does not saturate in small herds the way "ever free" does.

suppressPackageStartupMessages({library(SimInf); library(dplyr)})

stoch_scenario <- function(R0, N, e_s = 0.58, e_i = 0.74, D = Inf, p = 1, k = 20,
                           interval = NA, cv = 1, tau_days = 0, m_calf = 1,
                           reps = 200, years = 100, u = U_MORT,
                           detail = FALSE, tstep = 30, return_traj = FALSE) {
  Vn <- paste0("V", 1:k); comp <- c("S", "V0", Vn, "I", "IV")
  Ntot <- paste(comp, collapse = "+")
  inf <- sprintf("b0*(I + (1-e_i)*IV)/(%s)", Ntot)
  tr <- c("@ -> (1-p)*br*N0 -> S",
          "@ -> p*br*N0*has_onset -> V0",
          "@ -> p*br*N0*(1-has_onset) -> V1",
          sprintf("S -> S > 0 ? %s*S : 0 -> I", inf),
          sprintf("V0 -> V0 > 0 ? m_calf*%s*V0 : 0 -> I", inf),
          "V0 -> V0 > 0 ? g*V0 : 0 -> V1",
          "S -> S > 0 ? S*mu : 0 -> @", "V0 -> V0 > 0 ? V0*mu : 0 -> @",
          "I -> I > 0 ? I*mu : 0 -> @", "IV -> IV > 0 ? IV*mu : 0 -> @")
  for (j in 1:k) {
    nxt <- if (j < k) Vn[j + 1] else "S"
    tr <- c(tr, sprintf("%s -> %s > 0 ? (1-e_s)*%s*%s : 0 -> IV", Vn[j], Vn[j], inf, Vn[j]),
            sprintf("%s -> %s > 0 ? %s*mu : 0 -> @", Vn[j], Vn[j], Vn[j]),
            sprintf("%s -> %s > 0 ? %s*w : 0 -> %s", Vn[j], Vn[j], Vn[j], nxt))
  }
  mu <- u / 365
  R0v <- rep(R0, reps); Nv <- rep(N, reps)
  S0 <- round(Nv / R0v); I0 <- Nv - S0
  u0 <- as.data.frame(matrix(0L, reps, length(comp), dimnames = list(NULL, comp)))
  u0$S <- S0; u0$I <- I0
  ldata <- data.frame(N0 = Nv, b0 = R0v * mu, br = mu, mu = mu)
  gdata <- c(p = p, e_s = e_s, e_i = e_i, w = if (is.infinite(D)) 0 else k / (365 * D),
             has_onset = as.numeric(tau_days > 0), g = if (tau_days > 0) 1 / tau_days else 0, m_calf = m_calf)
  E <- matrix(0, length(comp), 1, dimnames = list(comp, "1")); Nm <- E
  iS <- match("S", comp); iV1 <- match("V1", comp)
  E[iS, 1] <- 1; Nm[iS, 1] <- iV1 - iS
  if (k > 1) for (j in 2:k) { i <- match(Vn[j], comp); E[i, 1] <- 1; Nm[i, 1] <- -(j - 1) }
  events <- NULL
  if (!is.na(interval)) {
    ev_days <- as.integer(round(seq(365 * interval, 365 * years, by = 365 * interval)))
    events <- data.frame(event = "intTrans", time = rep(ev_days, each = reps),
                         node = as.integer(rep(1:reps, times = length(ev_days))),
                         dest = 0L, n = 0L, proportion = cv, select = 1L, shift = 1L)
  }
  m <- mparse(transitions = tr, compartments = comp, ldata = ldata, gdata = gdata,
              u0 = u0, tspan = seq(1, 365 * years, by = tstep), events = events, E = E, N = Nm)
  x <- trajectory(SimInf::run(m))
  x$inf <- x$I + x$IV; x$tot <- rowSums(x[, comp]); x$yr <- x$time / 365
  if (return_traj) {
    return(x %>% filter(tot > 0) %>% group_by(yr) %>%
             summarise(mean_prev = mean(inf / tot), p_free = mean(inf == 0), .groups = "drop"))
  }
  # A small closed herd can die out entirely, which leaves prevalence
  # undefined. Record that explicitly and average prevalence over the
  # timepoints at which the herd still existed.
  s <- x %>% group_by(node) %>%
    summarise(T_free = if (any(inf == 0)) min(yr[inf == 0]) else Inf,
              extinct_20 = any(tot[yr <= 20] == 0),
              extinct_50 = any(tot[yr <= 50] == 0),
              mean_prev_20 = mean((inf / tot)[yr <= 20 & tot > 0]),
              mean_prev_50 = mean((inf / tot)[yr <= 50 & tot > 0]), .groups = "drop")
  if (detail) return(as.data.frame(s))
  s$T_free
}

# Summary used in the memo figures: median years to no infected animal
# (Inf if fewer than half the herds get there within the run) and the share
# free by year 50.
summ_T <- function(et) c(median_T = median(et), p50 = mean(et <= 50), p20 = mean(et <= 20))
