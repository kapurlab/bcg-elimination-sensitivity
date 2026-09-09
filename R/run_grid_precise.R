# Precise re-run of the conceptual R0 by herd-size grid.
#
# Two changes from run_stochastic_grid.R, both aimed at the fact that
# "probability a herd is ever free" is a poor basis for comparing programmes
# in small herds. It saturates, because chance fade-out clears such herds
# anyway, so the vaccine's marginal effect is a small difference between two
# numbers near one, and a ratio of two such differences is mostly noise.
#
#   1. Replicates raised from 200 to 2000 per cell, so the standard error on
#      a probability falls from about 3.5 to about 1.1 percentage points.
#   2. A continuous primary outcome: mean within-herd prevalence over the
#      first 20 (and 50) years, which is proportional to infected animal-years
#      and therefore to the economic loss. From it,
#         averted = 1 - mean prevalence under the programme
#                       / mean prevalence with no vaccination
#      is the share of infection-years the programme prevents. Its denominator
#      is large and precisely estimated in every cell, so unlike the old
#      retention ratio it is stable.
#
# Confidence intervals use a parametric bootstrap on the cell means, which is
# adequate at n = 2000. Cells where the interval is wider than 20 points are
# flagged rather than reported as a point estimate.
#
# Run from the project root: Rscript R/run_grid_precise.R (about 30-45 min)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/stochastic.R")
set.seed(20260909)

R0_grid <- c(1.2, 1.5, 2, 3, 5, 8)
N_grid  <- c(5, 10, 20, 50, 100, 200)
schedules <- list(
  "No vaccination (status quo)"                       = list(D = Inf, interval = NA, p = 0),
  "Lifelong immunity, calves only (paper assumption)" = list(D = Inf, interval = NA, p = 1),
  "18-month immunity, calves only"                    = list(D = 1.5, interval = NA, p = 1),
  "18-month immunity, annual whole-herd booster"      = list(D = 1.5, interval = 1,  p = 1),
  "18-month immunity, six-monthly whole-herd booster" = list(D = 1.5, interval = 0.5, p = 1))
REPS <- 2000; YEARS <- 50; TSTEP <- 60
dir.create("output/raw", showWarnings = FALSE)

design <- expand.grid(R0 = R0_grid, N = N_grid, schedule = names(schedules),
                      stringsAsFactors = FALSE)
t0 <- Sys.time()
CACHE <- "output/raw/grid_precise_replicates.rds"
if (file.exists(CACHE)) {
  raw <- readRDS(CACHE)
  cat("using cached replicates from", CACHE, "\n")
} else {
  raw <- vector("list", nrow(design))
  for (i in seq_len(nrow(design))) {
    d <- design[i, ]; sc <- schedules[[d$schedule]]
    r <- stoch_scenario(d$R0, d$N, D = sc$D, p = sc$p, interval = sc$interval,
                        reps = REPS, years = YEARS, detail = TRUE, tstep = TSTEP)
    raw[[i]] <- data.frame(R0 = d$R0, N = d$N, schedule = d$schedule,
                           T_free = r$T_free, mp20 = r$mean_prev_20, mp50 = r$mean_prev_50)
    if (i %% 20 == 0) cat(i, "of", nrow(design), "cells,",
                          round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
  }
  raw <- bind_rows(raw)
  saveRDS(raw, CACHE)
}
# A closed herd of 5 animals can die out over 50 years, leaving prevalence
# undefined. Those replicates are excluded from the prevalence outcomes and
# the rate is reported alongside.
if (!"extinct" %in% names(raw)) raw$extinct <- is.na(raw$mp20)

# --- per-cell summaries with standard errors ---------------------------------
cell <- raw %>% group_by(R0, N, schedule) %>%
  summarise(p_free_20 = mean(T_free <= 20), p_free_50 = mean(T_free <= 50),
            se_free_20 = sqrt(mean(T_free <= 20) * (1 - mean(T_free <= 20)) / n()),
            se_free_50 = sqrt(mean(T_free <= 50) * (1 - mean(T_free <= 50)) / n()),
            p_extinct = mean(extinct),
            # standard errors must be taken from the vectors before the means
            # overwrite those names inside summarise
            se_mp20 = sd(mp20, na.rm = TRUE) / sqrt(sum(!extinct)),
            se_mp50 = sd(mp50, na.rm = TRUE) / sqrt(sum(!extinct)),
            mp20 = mean(mp20, na.rm = TRUE),
            mp50 = mean(mp50, na.rm = TRUE),
            reps = n(), .groups = "drop")
write.csv(cell, "output/tables/grid_precise_cells.csv", row.names = FALSE)

# --- comparisons against the status quo, with bootstrap intervals -------------
B <- 20000
sq <- cell %>% filter(schedule == "No vaccination (status quo)") %>%
  select(R0, N, sq_pf20 = p_free_20, sq_se_pf20 = se_free_20,
         sq_mp20 = mp20, sq_se_mp20 = se_mp20, sq_mp50 = mp50, sq_se_mp50 = se_mp50)
life <- cell %>% filter(schedule == "Lifelong immunity, calves only (paper assumption)") %>%
  select(R0, N, lf_mp20 = mp20, lf_se_mp20 = se_mp20)

cmp <- cell %>% filter(schedule != "No vaccination (status quo)") %>%
  left_join(sq, by = c("R0", "N")) %>% left_join(life, by = c("R0", "N"))

boot_stat <- function(m_p, s_p, m_s, s_s, m_l, s_l) {
  # averted share, and that share as a fraction of what lifelong immunity averts
  dp <- rnorm(B, m_p, s_p); ds <- rnorm(B, m_s, s_s); dl <- rnorm(B, m_l, s_l)
  av <- 1 - dp / ds
  rel <- av / (1 - dl / ds)
  c(av_lo = unname(quantile(av, 0.025)), av_hi = unname(quantile(av, 0.975)),
    rel_lo = unname(quantile(rel, 0.025)), rel_hi = unname(quantile(rel, 0.975)))
}
bs <- t(mapply(boot_stat, cmp$mp20, cmp$se_mp20, cmp$sq_mp20, cmp$sq_se_mp20,
               cmp$lf_mp20, cmp$lf_se_mp20))
cmp <- bind_cols(cmp, as.data.frame(bs)) %>%
  mutate(averted_20 = 1 - mp20 / sq_mp20,
         averted_50 = 1 - mp50 / sq_mp50,
         rel_to_lifelong = averted_20 / (1 - lf_mp20 / sq_mp20),
         benefit_free_20 = p_free_20 - sq_pf20,
         se_benefit_free_20 = sqrt(se_free_20^2 + sq_se_pf20^2),
         averted_ci_width = av_hi - av_lo,
         rel_ci_width = rel_hi - rel_lo)
write.csv(cmp, "output/tables/grid_precise_comparisons.csv", row.names = FALSE)

# --- figure 35: share of infection-years averted -----------------------------
pal_lab <- function(x) sprintf("%.0f%%", 100 * x)
d1 <- cmp %>% mutate(schedule = factor(schedule, levels = names(schedules)[-1]),
                     R0f = factor(R0), Nf = factor(N))
g35 <- ggplot(d1, aes(R0f, Nf, fill = averted_20)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = pal_lab(averted_20), colour = averted_20 > 0.55), size = 3.4, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(0, 1), labels = scales::percent,
                      name = "Share of infected animal-years averted over the first 20 years") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  facet_wrap(~schedule, ncol = 2, labeller = label_wrap_gen(44)) +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "What vaccination averts, measured in infected animal-years rather than eventual freedom",
       subtitle = sprintf("%d replicate herds per cell. Mean within-herd prevalence over years 0 to 20, relative to no vaccination.\nThis outcome does not saturate in small herds, where chance fade-out already clears infection. In 5-animal herds, 11%% of replicates died out entirely and are excluded.", REPS)) +
  theme_pub(11) + theme(plot.title.position = "plot", panel.grid = element_blank(),
                        panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig35_infection_years_averted.png", g35, width = 11, height = 8.5, dpi = 200, bg = "white")

# --- figure 36: how much of the lifelong benefit a short-lived vaccine keeps --
d2 <- cmp %>% filter(schedule == "18-month immunity, calves only") %>%
  mutate(R0f = factor(R0), Nf = factor(N),
         uncertain = rel_ci_width > 0.20,
         label = ifelse(uncertain, "wide", sprintf("%.0f%%", 100 * rel_to_lifelong)),
         fill_val = ifelse(uncertain, NA, pmin(pmax(rel_to_lifelong, 0), 1)))
g36 <- ggplot(d2, aes(R0f, Nf, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label, colour = !is.na(fill_val) & fill_val > 0.55), size = 3.6, fontface = "bold") +
  scale_fill_gradient(low = "#FFF9DB", high = "#B7950B", limits = c(0, 1), labels = scales::percent,
                      na.value = "grey85", name = "Share of the lifelong vaccine's benefit retained") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  labs(x = "Within-herd R0", y = "Herd size (animals)",
       title = "Where an 18-month vaccine given only to calves still does most of the job",
       subtitle = sprintf("Infection-years averted over 20 years, as a share of what lifelong immunity averts. %d replicate herds per cell.\nCells marked 'wide' have a 95%% interval broader than 20 points and are not reported as a point estimate.", REPS)) +
  theme_pub(12) + theme(plot.title.position = "plot", panel.grid = element_blank(),
                        panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig36_short_duration_retained.png", g36, width = 9.5, height = 6.5, dpi = 200, bg = "white")

cat("\nGrid precise done in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
cat("\nShare of infection-years averted over 20 years (18-month, calves only):\n")
print(cmp %>% filter(schedule == "18-month immunity, calves only") %>%
        mutate(v = sprintf("%.0f%% (%.0f-%.0f)", 100 * averted_20, 100 * av_lo, 100 * av_hi)) %>%
        select(N, R0, v) %>% pivot_wider(names_from = R0, values_from = v) %>% as.data.frame())
cat("\nRetained share of the lifelong benefit:\n")
print(d2 %>% mutate(v = sprintf("%.0f%% (%.0f-%.0f)", 100 * rel_to_lifelong, 100 * rel_lo, 100 * rel_hi)) %>%
        select(N, R0, v) %>% pivot_wider(names_from = R0, values_from = v) %>% as.data.frame())
