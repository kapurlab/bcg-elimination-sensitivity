# Simplified replacement-policy figure and the archetype table for the memo,
# drawn from tables already produced by run_alt_views.R and run_herd.R.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/root.R"); source("R/model.R")
tr <- read.csv("output/tables/replacement_trajectories.csv")
keep_pol <- c("All replacements home-bred", "Half purchased, source 5%", "Half purchased, source 20%")
tp <- tr %>% filter(policy %in% keep_pol, NA_adults %in% c(20, 100)) %>%
  mutate(policy = factor(policy, levels = keep_pol),
         programme = factor(programme, levels = c("No vaccination", "18-month immunity, annual whole-herd booster")),
         size = factor(sprintf("%d adults (about %d head)", NA_adults, ifelse(NA_adults == 20, 44, 221)),
                       levels = c("20 adults (about 44 head)", "100 adults (about 221 head)")))
g <- ggplot(tp, aes(yr, prev_med, colour = programme, fill = programme)) +
  geom_ribbon(aes(ymin = prev_q25, ymax = prev_q75), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("grey40", "#1E8449"), name = NULL) + scale_fill_manual(values = c("grey40", "#1E8449"), name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(policy ~ size, labeller = label_wrap_gen(26)) +
  labs(x = "Years since the start of the programme", y = "Herd prevalence (median and interquartile band, 200 herds)",
       title = "Where replacements come from decides what vaccination can achieve",
       subtitle = "Age-structured herds at R0 3; adults leave at 0.2 per year and are replaced at once, home-bred or bought from a source of the stated prevalence.\nPurchased animals are unvaccinated. Vaccination lowers prevalence but cannot remove the floor that purchases set.") +
  theme_pub(11) + theme(plot.title.position = "plot")
ggsave("output/figs/fig34_memo_replacement.png", g, width = 11, height = 8, dpi = 200, bg = "white")

hr <- read.csv("output/tables/herd_replacement.csv") %>% filter(R0 == 3)
pick <- tribble(~archetype, ~NA_adults, ~policy,
  "Smallholder (5 adults), replacements home-bred", 5, "All replacements home-bred",
  "Smallholder (5 adults), half bought from a 20% source", 5, "Half purchased, source prevalence 20%",
  "Commercial (100 adults), replacements home-bred", 100, "All replacements home-bred",
  "Commercial (100 adults), half bought from a 5% source", 100, "Half purchased, source prevalence 5%",
  "Commercial (100 adults), half bought from a 20% source", 100, "Half purchased, source prevalence 20%")
vac <- hr %>% filter(programme != "No vaccination"); sq <- hr %>% filter(programme == "No vaccination")
tab <- pick %>% left_join(vac %>% select(NA_adults, policy, p_free_20, p_free_50, prev_50), by = c("NA_adults", "policy")) %>%
  left_join(sq %>% select(NA_adults, policy, sq_prev_50 = prev_50), by = c("NA_adults", "policy")) %>%
  transmute(Archetype = archetype, `Free by year 20` = sprintf("%.0f%%", 100 * p_free_20), `Free by year 50` = sprintf("%.0f%%", 100 * p_free_50),
            `Prevalence held at year 50` = sprintf("%.1f%%", 100 * prev_50), `Prevalence at year 50 without vaccination` = sprintf("%.0f%%", 100 * sq_prev_50))
write.csv(tab, "output/tables/memo_archetype_table.csv", row.names = FALSE)
cat("Memo extras done\n"); print(as.data.frame(tab))
