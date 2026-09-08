# Persistent non-responders against per-dose random non-response.
# Run from the project root: Rscript R/run_nonresponders.R

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/profile.R")
herds <- load_herds("data"); R0 <- median(herds$R0); N <- 100
VE <- 1 - (1 - BASE$e_s) * (1 - BASE$e_i)
K <- 40; SL <- 0.1
prof24 <- make_profile(K, SL, 30, 2.0, 2.05)      # near-step, 24 months
q_grid <- c(0, 0.1, 0.2, 0.3, 0.4)

run_q <- function(q, persistent, interval = 1) {
  s <- if (persistent) simulate_profile(R0, N, BASE$e_s, BASE$e_i, prof24, stage_len = SL, interval = interval, q = q)
       else simulate_profile(R0, N, BASE$e_s, BASE$e_i, prof24, stage_len = SL, interval = interval, p = 1 - q, cv = 1 - q)
  crossing_time(s$times, s$prev[, 1], 1e-3)
}
design <- expand.grid(q = q_grid, type = c("Persistent non-responders (never protected by any dose)",
                                            "Random non-response (each dose fails independently)"),
                      stringsAsFactors = FALSE)
design$T <- mapply(function(q, ty) run_q(q, grepl("Persistent", ty)), design$q, design$type)
write.csv(design, "output/tables/nonresponders.csv", row.names = FALSE)

# analytic ceiling on persistent non-response, lifelong immunity, full coverage
R0s <- exp(seq(log(1.2), log(18), length.out = 200))
ceiling <- data.frame(R0 = R0s, q_max = pmax(0, 1 - (1 - 1 / R0s) / VE))
herd_pts <- data.frame(R0 = herds$R0, q_max = pmax(0, 1 - (1 - 1 / herds$R0) / VE))
qR0 <- quantile(herds$R0, c(0.25, 0.5, 0.75, 0.9))

gA <- ggplot(ceiling, aes(R0, q_max)) +
  geom_line(linewidth = 1.1, colour = "#7B1E3A") +
  geom_point(data = herd_pts, colour = "#7FB3D5", size = 2, alpha = 0.8) +
  geom_vline(xintercept = qR0, linetype = 3, colour = "grey55") +
  annotate("text", x = qR0, y = 0.62, label = c("25th", "median", "75th", "90th"), size = 3, colour = "grey30", vjust = 0) +
  scale_x_log10(breaks = c(1.5, 2, 3, 5, 8, 12, 17)) + scale_y_continuous(labels = scales::percent, limits = c(0, 0.66)) +
  labs(x = "Herd R0 (log scale)", y = "Largest tolerable share of persistent non-responders",
       title = "A. Ceiling on persistent non-responders for elimination to remain possible",
       subtitle = "Lifelong immunity, every other animal protected, 89% total efficacy. Points: the 57 study herds.") +
  theme_pub(12) + theme(plot.title.position = "plot")
ggsave("output/figs/fig19a_nonresponder_ceiling.png", gA, width = 9, height = 5.5, dpi = 200, bg = "white")

d <- design %>% mutate(type = factor(type, levels = rev(unique(type))), qlab = sprintf("%d%%", round(100 * q)),
                       never = !is.finite(T), fill_val = ifelse(never, NA, pmin(T, 120)),
                       label = ifelse(never, "never", sprintf("%.0f", T)))
gB <- ggplot(d, aes(qlab, type, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = label, colour = never | fill_val > 80), size = 4.6, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(25, 120), na.value = "grey88",
                      name = "Years to herd prevalence < 0.1%", breaks = c(25, 50, 75, 100, 120), labels = c("25", "50", "75", "100", "120+")) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, 28)) +
  labs(x = "Share of animals not responding", y = NULL,
       title = "B. Non-response under annual whole-herd revaccination, 24-month immunity, median herd",
       subtitle = "Random non-response is re-drawn at every dose; persistent non-responders are never protected.") +
  theme_pub(12) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig19b_nonresponder_tiles.png", gB, width = 10, height = 4.5, dpi = 200, bg = "white")
cat("Non-responders done\n"); print(design %>% mutate(T = round(T)) %>% pivot_wider(names_from = q, values_from = T))
cat("ceiling at herd quantiles:\n"); print(round(pmax(0, 1 - (1 - 1 / qR0) / VE), 2))
