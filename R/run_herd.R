# Age-structured stochastic herd with replacement policy: home-bred against
# purchased replacements from a source population of given prevalence, across
# adult herd size and within-herd R0, with and without vaccination.
# Run from the project root: Rscript R/run_herd.R (about 10 minutes)

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R"); source("R/herd.R")
set.seed(20260908)
mod <- build_herd_model(k = 10)
NA_grid <- c(5, 20, 100); R0_grid <- c(1.5, 3, 5)
policies <- tribble(
  ~policy, ~h, ~pi_src,
  "All replacements home-bred", 1, 0,
  "Half purchased, source prevalence 5%", 0.5, 0.05,
  "Half purchased, source prevalence 20%", 0.5, 0.20,
  "Half purchased, source prevalence 50%", 0.5, 0.50,
  "All purchased, source prevalence 5%", 0, 0.05,
  "All purchased, source prevalence 20%", 0, 0.20,
  "All purchased, source prevalence 50%", 0, 0.50)
programmes <- tribble(
  ~programme, ~p, ~interval,
  "No vaccination", 0, NA,
  "18-month immunity, calves at birth plus annual whole-herd booster", 1, 1)
REPS <- 200; YEARS <- 60

design <- crossing(NA_adults = NA_grid, R0 = R0_grid, policies, programmes)
res <- bind_rows(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, ]
  r <- run_herd(mod, R0 = d$R0, NA_adults = d$NA_adults, p = d$p, interval = d$interval,
                h = d$h, pi_src = d$pi_src, reps = REPS, years = YEARS)
  data.frame(d, p_free_20 = mean(r$free_20), p_free_50 = mean(r$free_50),
             prev_20 = mean(r$prev_20), prev_50 = mean(r$prev_50), median_T = median(r$T_free),
             herd_size_50 = mean(r$N_50))
}))
write.csv(res, "output/tables/herd_replacement.csv", row.names = FALSE)

# --- figures -------------------------------------------------------------------
sizes <- res %>% filter(programme != "No vaccination", h == 1) %>% group_by(NA_adults) %>%
  summarise(lab = sprintf("%d adults (about %.0f head)", first(NA_adults), mean(herd_size_50)), .groups = "drop") %>% arrange(NA_adults)
d <- res %>% mutate(policy = factor(policy, levels = rev(policies$policy)), R0f = factor(R0),
                    programme = factor(programme, levels = programmes$programme),
                    size_lab = factor(sizes$lab[match(NA_adults, sizes$NA_adults)], levels = sizes$lab))

g24 <- ggplot(d, aes(R0f, policy, fill = p_free_50)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * p_free_50), colour = p_free_50 > 0.55), size = 3.2, fontface = "bold") +
  scale_fill_gradient(low = "#FBEAEA", high = "#7B1E3A", limits = c(0, 1), labels = scales::percent,
                      name = "Probability the herd has no infected animal at year 50") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, 26)) +
  facet_grid(programme ~ size_lab, labeller = label_wrap_gen(28)) +
  labs(x = "Within-herd R0", y = NULL,
       title = "Replacement policy decides what vaccination can achieve: probability a herd is free of infection at year 50",
       subtitle = sprintf("Age-structured stochastic herds, %d replicates per cell, 20-year burn-in at the endemic level. Adults leave at 0.2 per year and are replaced at once, home-bred or purchased.\nPurchased animals are unvaccinated and infected with the source prevalence. Direct efficacy 58%%, indirect 74%%.", REPS)) +
  theme_pub(10) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig24_replacement_policy_free50.png", g24, width = 13, height = 9, dpi = 200, bg = "white")

g25 <- ggplot(d, aes(R0f, policy, fill = pmin(prev_50, 0.5))) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(prev_50 < 0.0005, "0", sprintf("%.1f%%", 100 * prev_50)), colour = prev_50 > 0.28), size = 3.2, fontface = "bold") +
  scale_fill_gradient(low = "#FFF9DB", high = "#B7950B", limits = c(0, 0.5), labels = scales::percent,
                      name = "Mean herd prevalence at year 50 (capped at 50%)") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, 26)) +
  facet_grid(programme ~ size_lab, labeller = label_wrap_gen(28)) +
  labs(x = "Within-herd R0", y = NULL,
       title = "Prevalence held at year 50 under each replacement policy",
       subtitle = "Same runs as the previous figure. Purchased infected animals set the floor that vaccination cannot remove.") +
  theme_pub(10) + theme(plot.title.position = "plot", panel.grid = element_blank(), panel.border = element_blank(), legend.key.width = unit(1.6, "cm"))
ggsave("output/figs/fig25_replacement_policy_prev50.png", g25, width = 13, height = 9, dpi = 200, bg = "white")
cat("Herd replacement done\n")
print(res %>% filter(programme != "No vaccination") %>% mutate(v = sprintf("%.0f%% / %.1f%%", 100 * p_free_50, 100 * prev_50)) %>%
        select(NA_adults, R0, policy, v) %>% pivot_wider(names_from = R0, values_from = v), n = 40, width = 200)
