# One memo figure: years to elimination in the median herd for named
# scenarios, one lever at a time. Run from the project root:
# Rscript R/run_memo_figure.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("R/root.R")
source("R/model.R")
source("R/revaccination.R")
source("R/onset.R")
herds <- load_herds("data")
R0 <- median(herds$R0)
N <- 100
HOR <- 200

T_base <- function(e_s = BASE$e_s, e_i = BASE$e_i, D = Inf, p = 1)
  T_single(R0, e_s, e_i, D, p, N = N, horizon = HOR)
T_rev <- function(D, k = 20, revacc = TRUE, redose_V = FALSE) {
  s <- simulate_revacc(R0, N, BASE$e_s, BASE$e_i, D = D, k = k, revacc = revacc,
                       redose_V = redose_V, horizon = HOR)
  crossing_time(s$times, s$prev[, 1], 1e-3)
}
T_on <- function(tau_days, m) {
  s <- simulate_onset(R0, N, BASE$e_s, BASE$e_i, D = Inf, tau_days = tau_days, m = m,
                      k = 1, horizon = HOR)
  crossing_time(s$times, s$prev[, 1], 1e-3)
}

sc <- tribble(
  ~lever, ~scenario, ~T,
  "Paper base case", "89% total efficacy, lifelong protection, all calves", T_base(),
  "Peak efficacy", "95% total efficacy", T_base(0.75, 0.80),
  "Peak efficacy", "100% total efficacy", T_base(0.95, 0.95),
  "Coverage of calves", "90% of calves vaccinated", T_base(p = 0.9),
  "Coverage of calves", "80% of calves vaccinated", T_base(p = 0.8),
  "Coverage of calves", "70% of calves vaccinated", T_base(p = 0.7),
  "Duration, calves only", "10 years of protection", T_rev(10, revacc = FALSE),
  "Duration, calves only", "5 years of protection", T_rev(5, revacc = FALSE),
  "Duration, calves only", "18 months of protection", T_rev(1.5, revacc = FALSE),
  "Duration, annual boost of lapsed animals", "36 months of protection", T_rev(3),
  "Duration, annual boost of lapsed animals", "24 months of protection", T_rev(2),
  "Duration, annual boost of lapsed animals", "18 months of protection", T_rev(1.5),
  "Duration, annual boost of lapsed animals", "12 months of protection", T_rev(1),
  "Duration, annual boost of all uninfected animals", "18 months of protection", T_rev(1.5, redose_V = TRUE),
  "Duration, annual boost of all uninfected animals", "12 months of protection", T_rev(1, redose_V = TRUE),
  "Duration, annual boost of all uninfected animals", "6 months of protection", T_rev(0.5, redose_V = TRUE),
  "30-day window before protection", "calf exposure equal to adults", T_on(30, 1),
  "30-day window before protection", "calf exposure 3 x adults", T_on(30, 3),
  "30-day window before protection", "calf exposure 10 x adults", T_on(30, 10))
write.csv(sc, "output/tables/memo_scenarios.csv", row.names = FALSE)

lev <- unique(sc$lever)
pal <- c("Paper base case" = "grey30", "Peak efficacy" = "#7FB3D5",
         "Coverage of calves" = "#F9E79F", "Duration, calves only" = "#F1948A",
         "Duration, annual boost of lapsed animals" = "#A9DFBF",
         "Duration, annual boost of all uninfected animals" = "#1E8449",
         "30-day window before protection" = "#BB8FCE")
CAP <- 150
d <- sc %>% mutate(lever = factor(lever, levels = lev),
                   never = !is.finite(T) | T > CAP,
                   x = ifelse(never, CAP, T),
                   label = ifelse(never, "never", sprintf("%.0f y", T)),
                   row = factor(paste(lever, scenario), levels = rev(paste(sc$lever, sc$scenario))))
base <- sc$T[1]
g <- ggplot(d, aes(x = x, y = row)) +
  geom_vline(xintercept = base, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = 0, xend = x, yend = row, colour = lever, alpha = never), linewidth = 2.2, lineend = "round") +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4), guide = "none") +
  geom_point(data = filter(d, never), shape = 4, size = 3.5, stroke = 1.2, colour = "grey20") +
  geom_text(aes(label = label), hjust = -0.25, size = 3.4, colour = "grey20") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_y_discrete(labels = setNames(d$scenario, d$row)) +
  scale_x_continuous(limits = c(0, CAP + 25), breaks = c(0, 25, 50, 75, 100, 125, CAP),
                     labels = c("0", "25", "50", "75", "100", "125", "no\nelimination"),
                     expand = expansion(mult = c(0, 0))) +
  facet_grid(lever ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = label_wrap_gen(22)) +
  labs(x = "Years from the start of vaccination to herd prevalence below 0.1%", y = NULL,
       title = "What each lever does to time to elimination in a typical Ethiopian dairy herd",
       subtitle = sprintf("Median herd, R0 = %.1f. Each row changes one thing from the paper's base case (dashed line, %.0f years).\nCrosses: the reproduction number under vaccination stays above 1, so prevalence never reaches the threshold.", R0, base)) +
  theme_pub(11) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        strip.background = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9.5), panel.spacing.y = unit(4, "pt"))
ggsave("output/figs/fig14_memo_scenarios.png", g, width = 11, height = 8.5, dpi = 200, bg = "white")
cat("Memo figure done\n"); print(as.data.frame(sc))
