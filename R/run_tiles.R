# Tile figure: booster schedule by duration of protection, median herd.
# Panels 1 to 3: booster duration equal to, twice, or independent of the
# birth-dose duration. Panel 4: coverage of the booster campaigns.
# Run from the project root: Rscript R/run_tiles.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("R/root.R")
source("R/model.R")
source("R/boost.R")
herds <- load_herds("data")
R0 <- median(herds$R0)
N <- 100

T_boost <- function(D1, D2, interval, cv = 1) {
  s <- simulate_boost(R0, N, BASE$e_s, BASE$e_i, D1 = D1, D2 = D2, interval = interval, cv = cv)
  crossing_time(s$times, s$prev[, 1], 1e-3)
}

sched <- data.frame(schedule = c("No booster", "Booster every 24 months", "Booster annually", "Booster every 6 months"),
                    interval = c(NA, 2, 1, 0.5))
D1s <- c(0.5, 1, 1.5, 2)
D1lab <- function(D) sprintf("%d months", round(D * 12))

grid_panels <- bind_rows(
  merge(sched, data.frame(D1 = D1s)) %>% mutate(panel = "A. Booster lasts as long as the birth dose", D2 = D1),
  merge(sched, data.frame(D1 = D1s)) %>% mutate(panel = "B. Booster lasts twice as long as the birth dose", D2 = 2 * D1),
  merge(sched, data.frame(D1 = D1s)) %>% mutate(panel = "C. Booster lasts 5 years whatever the birth dose", D2 = 5)) %>%
  mutate(cv = 1, column = D1lab(D1))
cov_panel <- merge(sched %>% filter(!is.na(interval)), data.frame(cv = c(0.5, 0.7, 0.9, 1))) %>%
  mutate(panel = "D. Campaign coverage, birth dose and booster both 12 months", D1 = 1, D2 = 1,
         column = sprintf("%d%% of herd", round(100 * cv)))
# single grid: one duration of immunity whichever dose produced it
single <- merge(sched, data.frame(D1 = c(0.5, 1, 2, 10))) %>%
  mutate(panel = "Single grid", D2 = D1, cv = 1,
         column = ifelse(D1 >= 5, sprintf("%d years", D1), D1lab(D1)))
design <- bind_rows(grid_panels, cov_panel, single)

design$T <- mapply(function(D1, D2, iv, cv) T_boost(D1, D2, iv, cv), design$D1, design$D2, design$interval, design$cv)
write.csv(design, "output/tables/tiles_booster_schedule.csv", row.names = FALSE)

# --- figure ------------------------------------------------------------------
d_all <- design %>%
  mutate(schedule = factor(schedule, levels = rev(sched$schedule)),
         column = factor(column, levels = unique(column)),
         never = !is.finite(T),
         fill_val = ifelse(never, NA, pmin(T, 120)),
         label = ifelse(never, "never", sprintf("%.0f", T)),
         panel = factor(panel, levels = unique(panel)))
d <- d_all %>% filter(panel != "Single grid")
base_T <- T_single(R0, BASE$e_s, BASE$e_i, Inf, 1, N = N)
g <- ggplot(d, aes(column, schedule, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = label, colour = never | fill_val > 80), size = 4.2, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(25, 120), na.value = "grey88",
                      name = "Years to herd prevalence < 0.1%", breaks = c(25, 50, 75, 100, 120),
                      labels = c("25", "50", "75", "100", "120+")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  facet_wrap(~panel, scales = "free_x", ncol = 2) +
  labs(x = NULL, y = NULL,
       title = "Booster schedule against duration of protection, median Ethiopian dairy herd",
       subtitle = sprintf("Every calf dosed at birth. Columns in A to C: duration of birth-dose protection. Grey: elimination impossible. Paper's lifelong-protection case: %.0f years.\nBoosters go to the whole herd, protected or not, and dosing a still-protected animal is assumed neutral. R0 = %.1f, direct efficacy 58%%, indirect 74%%.", base_T, R0)) +
  theme_pub(12) +
  theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(),
        legend.key.width = unit(1.6, "cm"), strip.text = element_text(hjust = 0),
        axis.text.x = element_text(size = 10))
ggsave("output/figs/fig15_booster_tiles.png", g, width = 12, height = 8, dpi = 200, bg = "white")

# --- single grid figure ---------------------------------------------------
ds <- d_all %>% filter(panel == "Single grid") %>%
  mutate(column = factor(column, levels = c("6 months", "12 months", "24 months", "10 years")),
         schedule = factor(as.character(schedule), levels = rev(c("No booster", "Booster every 6 months", "Booster annually", "Booster every 24 months"))))
g2 <- ggplot(ds, aes(column, schedule, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = never | fill_val > 80), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(25, 120), na.value = "grey88",
                      name = "Years to herd prevalence < 0.1%", breaks = c(25, 50, 75, 100, 120),
                      labels = c("25", "50", "75", "100", "120+")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  labs(x = "Duration of immunity after any dose", y = "Revaccination interval",
       title = "Revaccination interval against duration of immunity, median Ethiopian dairy herd",
       subtitle = sprintf("Every calf dosed at birth; boosters go to the whole herd. Grey: elimination impossible. Paper's lifelong case: %.0f years.\nR0 = %.1f, direct efficacy 58%%, indirect 74%%, herd turnover 0.27 per year.", base_T, R0)) +
  theme_pub(13) +
  theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(),
        legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig16_interval_by_duration.png", g2, width = 9.5, height = 6, dpi = 200, bg = "white")
cat("Tiles done\n")
print(d %>% select(panel, schedule, column, label) %>% tidyr::pivot_wider(names_from = column, values_from = label), n = 30, width = 200)
