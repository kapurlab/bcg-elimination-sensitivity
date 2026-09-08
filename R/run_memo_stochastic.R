# Stochastic versions of the two memo figures for one named herd: 100 head
# at the survey-median R0 (2.8), 200 replicate herds per scenario, outcome
# the median years until the herd has no infected animal.
# Run from the project root: Rscript R/run_memo_stochastic.R (about 5 min)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/stochastic.R")
set.seed(20260909)
herds <- load_herds("data"); R0 <- median(herds$R0); N <- 100; REPS <- 200; YEARS <- 100
run_scn <- function(...) summ_T(stoch_scenario(R0, N, reps = REPS, years = YEARS, ...))

# --- ladder ----------------------------------------------------------------------
sc <- tribble(
  ~lever, ~scenario, ~args,
  "Status quo", "No vaccination", list(p = 0),
  "Paper base case", "89% total efficacy, lifelong immunity, all calves", list(),
  "Peak efficacy, lifelong immunity", "95% total efficacy", list(e_s = 0.75, e_i = 0.80),
  "Peak efficacy, lifelong immunity", "100% total efficacy", list(e_s = 0.95, e_i = 0.95),
  "Coverage of calves at birth", "90% of calves vaccinated", list(p = 0.9),
  "Coverage of calves at birth", "80% of calves vaccinated", list(p = 0.8),
  "Coverage of calves at birth", "70% of calves vaccinated", list(p = 0.7),
  "Duration of immunity, birth dose only", "10 years of immunity", list(D = 10),
  "Duration of immunity, birth dose only", "5 years of immunity", list(D = 5),
  "Duration of immunity, birth dose only", "18 months of immunity", list(D = 1.5),
  "Duration of immunity, annual revaccination of the whole herd", "18 months of immunity", list(D = 1.5, interval = 1),
  "Duration of immunity, annual revaccination of the whole herd", "12 months of immunity", list(D = 1, interval = 1),
  "Duration of immunity, annual revaccination of the whole herd", "6 months of immunity", list(D = 0.5, interval = 1),
  "30-day window between birth and immunity, lifelong immunity", "calf exposure equal to adults", list(tau_days = 30, m_calf = 1),
  "30-day window between birth and immunity, lifelong immunity", "calf exposure 3 x adults", list(tau_days = 30, m_calf = 3),
  "30-day window between birth and immunity, lifelong immunity", "calf exposure 10 x adults", list(tau_days = 30, m_calf = 10))
out <- t(sapply(sc$args, function(a) do.call(run_scn, a)))
lad <- bind_cols(sc %>% select(lever, scenario), as.data.frame(out))
write.csv(lad, "output/tables/memo_scenarios_stochastic.csv", row.names = FALSE)

lev <- unique(lad$lever)
pal <- c("Status quo" = "grey60", "Paper base case" = "grey30", "Peak efficacy, lifelong immunity" = "#7FB3D5",
         "Coverage of calves at birth" = "#F9E79F", "Duration of immunity, birth dose only" = "#F1948A",
         "Duration of immunity, annual revaccination of the whole herd" = "#1E8449",
         "30-day window between birth and immunity, lifelong immunity" = "#BB8FCE")
CAP <- 100
d <- lad %>% mutate(lever = factor(lever, levels = lev), never = !is.finite(median_T),
                    x = ifelse(never, CAP, median_T),
                    label = ifelse(never, sprintf("fewer than half free by year 100 (%.0f%% by 50)", 100 * p50),
                                   sprintf("%.0f y (%.0f%% free by 50)", median_T, 100 * p50)),
                    row = factor(paste(lever, scenario), levels = rev(paste(lad$lever, lad$scenario))))
base <- lad$median_T[lad$lever == "Paper base case"]
g <- ggplot(d, aes(x = x, y = row)) +
  geom_vline(xintercept = base, linetype = 2, colour = "grey50") +
  geom_segment(aes(x = 0, xend = x, yend = row, colour = lever, alpha = never), linewidth = 2.2, lineend = "round") +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.4), guide = "none") +
  geom_point(data = filter(d, never), shape = 4, size = 3.5, stroke = 1.2, colour = "grey20") +
  geom_text(aes(label = label), hjust = -0.1, size = 3.2, colour = "grey20") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_y_discrete(labels = setNames(d$scenario, d$row)) +
  scale_x_continuous(limits = c(0, CAP + 60), breaks = c(0, 25, 50, 75, CAP), labels = c("0", "25", "50", "75", "100+"), expand = expansion(mult = 0)) +
  facet_grid(lever ~ ., scales = "free_y", space = "free_y", switch = "y", labeller = label_wrap_gen(26)) +
  labs(x = "Median years until the herd has no infected animal", y = NULL,
       title = "What each lever does to time to elimination in a 100-head herd at the survey-median R0",
       subtitle = sprintf("Stochastic model, %d replicate herds per row, R0 = %.1f, turnover 0.27 per year. Each row changes one thing from the paper's base case (dashed line, %.0f years).\nCrosses: fewer than half of the herds are free by year 100.", REPS, R0, base)) +
  theme_pub(11) +
  theme(plot.title.position = "plot", strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        strip.background = element_blank(), panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 9.5), panel.spacing.y = unit(4, "pt"))
ggsave("output/figs/fig32_memo_ladder_stochastic.png", g, width = 12, height = 9, dpi = 200, bg = "white")

# --- interval by duration grid ----------------------------------------------------
sched <- data.frame(schedule = c("No booster", "Booster every 6 months", "Booster annually", "Booster every 24 months"), interval = c(NA, 0.5, 1, 2))
dur <- data.frame(D = c(0.5, 1, 2, 10), dlab = c("6 months", "12 months", "24 months", "10 years"))
gd <- merge(sched, dur)
gres <- t(mapply(function(iv, D) run_scn(D = D, interval = iv), gd$interval, gd$D))
grid <- bind_cols(gd, as.data.frame(gres))
write.csv(grid, "output/tables/interval_duration_stochastic.csv", row.names = FALSE)
gt <- grid %>% mutate(schedule = factor(schedule, levels = rev(sched$schedule)), dlab = factor(dlab, levels = dur$dlab),
                      never = !is.finite(median_T), fill_val = ifelse(never, NA, pmin(median_T, 100)),
                      label = ifelse(never, sprintf("< half free\n(%.0f%% by 50)", 100 * p50), sprintf("%.0f\n(%.0f%% by 50)", median_T, 100 * p50)))
g2 <- ggplot(gt, aes(dlab, schedule, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = never | fill_val > 70), size = 4, fontface = "bold", lineheight = 0.9) +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(15, 100), na.value = "grey88",
                      name = "Median years until the herd has no infected animal", breaks = c(25, 50, 75, 100), labels = c("25", "50", "75", "100+")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  labs(x = "Duration of immunity after any dose", y = "Revaccination interval",
       title = "Revaccination interval against duration of immunity, 100-head herd at the survey-median R0",
       subtitle = sprintf("Stochastic model, %d replicate herds per cell, R0 = %.1f. Every calf dosed at birth; boosters go to the whole herd.\nGrey: fewer than half of herds free by year 100.", REPS, R0)) +
  theme_pub(13) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig33_interval_by_duration_stochastic.png", g2, width = 10, height = 6.2, dpi = 200, bg = "white")
cat("Memo stochastic done\n"); print(as.data.frame(lad)); print(as.data.frame(grid))
