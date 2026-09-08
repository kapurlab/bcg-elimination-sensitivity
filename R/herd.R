# Age-structured stochastic herd with a replacement policy (SimInf).
#
# Age classes: C calves (0 to 1 year), Y young stock (1 to 2 years), A adults.
# Within each class: S susceptible, I infected, IV infected while protected,
# V1..Vk protection stages (near-fixed duration for large k).
#
# Herd management:
#   births   : retained female calves at rate b per adult per year, into C.
#              A fraction p of them is dosed at birth (enter V1 of C).
#   ageing   : C -> Y at rate 1 per year; calf mortality mC.
#   adults   : fixed target number NA. Adults leave at rate cull per year.
#              Each leaver is replaced at once: with probability h by
#              promoting a young-stock animal (drawn from the Y compartments
#              in proportion), otherwise by a purchased animal that is
#              infected with probability pi_src and unvaccinated.
#   young    : young stock not promoted are sold at rate sig per year.
#   campaign : whole-herd booster every `interval` years moving S and all V
#              stages (all ages) to V1 of their class; dosing a protected
#              animal is neutral.
# Transmission is frequency dependent across the whole herd, with calves
# exposed at m_calf times the adult rate. b0 = R0 * cull, so R0 is the
# within-herd reproduction number for an adult infection lasting one
# adult residence (1 / cull years).

suppressPackageStartupMessages({library(SimInf); library(dplyr)})

build_herd_model <- function(k = 10) {
  ages <- c("C", "Y", "A")
  base <- c("S", "I", "IV", paste0("V", 1:k))
  comp <- as.vector(outer(base, ages, function(b, a) paste0(b, "_", a)))
  Vs <- function(a) paste0("V", 1:k, "_", a)
  tot <- function(a) paste(c(paste0("S_", a), paste0("I_", a), paste0("IV_", a), Vs(a)), collapse = "+")
  Itot <- paste(paste0("I_", ages), collapse = "+")
  IVtot <- paste(paste0("IV_", ages), collapse = "+")
  Ntot <- paste(sapply(ages, tot), collapse = "+")
  Atot <- tot("A"); Ytot <- tot("Y")
  lam <- sprintf("b0*(%s + (1-e_i)*(%s))/(%s)", Itot, IVtot, Ntot)
  mult <- c(C = "m_calf", Y = "1", A = "1")
  tr <- c(
    sprintf("@ -> b*(%s)*(1-p) -> S_C", Atot),
    sprintf("@ -> b*(%s)*p -> V1_C", Atot),
    # purchases replace the share (1 - h) of leavers, and also any leaver whose
    # home-bred replacement is unavailable because no young stock is present
    sprintf("@ -> ((1-h) + h*((%s) > 0 ? 0 : 1))*cull*NAd*pi_src -> I_A", Ytot),
    sprintf("@ -> ((1-h) + h*((%s) > 0 ? 0 : 1))*cull*NAd*(1-pi_src) -> S_A", Ytot))
  for (a in ages) {
    m <- mult[[a]]
    tr <- c(tr,
            sprintf("S_%s -> S_%s > 0 ? %s*%s*S_%s : 0 -> I_%s", a, a, m, lam, a, a))
    for (j in 1:k) {
      v <- sprintf("V%d_%s", j, a); nxt <- if (j < k) sprintf("V%d_%s", j + 1, a) else sprintf("S_%s", a)
      tr <- c(tr,
              sprintf("%s -> %s > 0 ? (1-e_s)*%s*%s*%s : 0 -> IV_%s", v, v, m, lam, v, a),
              sprintf("%s -> %s > 0 ? w*%s : 0 -> %s", v, v, v, nxt))
    }
    for (x in base) {
      cx <- sprintf("%s_%s", x, a)
      if (a == "C") tr <- c(tr, sprintf("%s -> %s > 0 ? %s*aC : 0 -> %s_Y", cx, cx, cx, x),
                            sprintf("%s -> %s > 0 ? %s*mC : 0 -> @", cx, cx, cx))
      if (a == "Y") tr <- c(tr, sprintf("%s -> (%s) > 0 && %s > 0 ? h*cull*NAd*%s/(%s) : 0 -> %s_A", cx, Ytot, cx, cx, Ytot, x),
                            sprintf("%s -> %s > 0 ? %s*sig : 0 -> @", cx, cx, cx))
      if (a == "A") tr <- c(tr, sprintf("%s -> %s > 0 ? %s*cull : 0 -> @", cx, cx, cx))
    }
  }
  # campaign: select column 1 picks S and V2..Vk of every age; shift moves them to V1 of that age
  E <- matrix(0, length(comp), 1, dimnames = list(comp, "1"))
  Nm <- matrix(0, length(comp), 1, dimnames = list(comp, "1"))
  for (a in ages) {
    iS <- match(sprintf("S_%s", a), comp); iV1 <- match(sprintf("V1_%s", a), comp)
    E[iS, 1] <- 1; Nm[iS, 1] <- iV1 - iS
    if (k > 1) for (j in 2:k) { i <- match(sprintf("V%d_%s", j, a), comp); E[i, 1] <- 1; Nm[i, 1] <- -(j - 1) }
  }
  list(comp = comp, transitions = tr, E = E, N = Nm, k = k)
}

# Run one scenario: burn-in without vaccination, then the programme.
# Returns per-replicate: first year free of infection (Inf if never),
# whether free at years 20 and 50, and prevalence at years 20 and 50.
run_herd <- function(model, R0, NA_adults, e_s = 0.58, e_i = 0.74, D = 1.5, p = 1,
                     interval = NA, h = 1, pi_src = 0, b = 0.45, cull = 0.2, sig = 0.3,
                     mC = 0.05, m_calf = 1, reps = 200, years = 100, burnin = 20) {
  k <- model$k; comp <- model$comp
  day <- 365
  gd <- function(pv, iv) c(b = b / day, p = pv, h = h, cull = cull / day, pi_src = pi_src,
                            e_s = e_s, e_i = e_i, w = if (is.infinite(D)) 0 else k / (D * day),
                            aC = 1 / day, mC = mC / day, sig = sig / day, m_calf = m_calf)
  ld <- data.frame(NAd = rep(NA_adults, reps), b0 = rep(R0 * cull / day, reps))
  u0 <- as.data.frame(matrix(0L, reps, length(comp), dimnames = list(NULL, comp)))
  I0 <- round(NA_adults * (1 - 1 / R0)); u0$I_A <- I0; u0$S_A <- NA_adults - I0
  u0$S_Y <- round(0.3 * NA_adults); u0$S_C <- round(0.4 * NA_adults)
  # burn-in: no vaccination, no campaigns
  m1 <- mparse(transitions = model$transitions, compartments = comp, ldata = ld, gdata = gd(0, NA),
               u0 = u0, tspan = c(1, burnin * day))
  x1 <- trajectory(run(m1)); x1 <- x1[x1$time == burnin * day, ]
  u1 <- x1[, comp]
  # programme
  tspan <- seq(1, years * day, by = 30)
  ev <- NULL
  if (!is.na(interval)) {
    days <- as.integer(round(seq(interval * day, years * day, by = interval * day)))
    ev <- data.frame(event = "intTrans", time = rep(days, each = reps), node = as.integer(rep(1:reps, times = length(days))),
                     dest = 0L, n = 0L, proportion = 1, select = 1L, shift = 1L)
  }
  m2 <- mparse(transitions = model$transitions, compartments = comp, ldata = ld, gdata = gd(p, interval),
               u0 = u1, tspan = tspan, events = ev, E = model$E, N = model$N)
  x <- trajectory(run(m2))
  inf_cols <- grep("^(I|IV)_", comp, value = TRUE)
  x$inf <- rowSums(x[, inf_cols]); x$N <- rowSums(x[, comp])
  x$yr <- x$time / day
  x %>% group_by(node) %>% summarise(
    T_free = if (any(inf == 0)) min(yr[inf == 0]) else Inf,
    free_20 = inf[which.min(abs(yr - 20))] == 0, free_50 = inf[which.min(abs(yr - 50))] == 0,
    prev_20 = inf[which.min(abs(yr - 20))] / pmax(1, N[which.min(abs(yr - 20))]),
    prev_50 = inf[which.min(abs(yr - 50))] / pmax(1, N[which.min(abs(yr - 50))]),
    N_50 = N[which.min(abs(yr - 50))], .groups = "drop")
}
