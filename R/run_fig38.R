# Figure 38: the 20 demographic strata of the benefit-timing analysis, drawn
# separately by herd size and within-herd R0.
#
# Plotting layer only. The model output is read from the tables written by
# R/run_benefit_timing.R and is not recomputed here, so this script can be
# iterated on freely without touching the simulation.
#
# Run from the project root: Rscript R/run_fig38.R (seconds)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(cowplot)
  library(ggrepel); library(RColorBrewer)
})
source("R/root.R"); source("R/model.R")

N_levels <- c(10, 20, 50, 100, 200)
tr <- read.csv("output/tables/benefit_timing.csv") %>% mutate(Nf = factor(N, levels = N_levels))
sm <- read.csv("output/tables/benefit_timing_summary.csv") %>% mutate(Nf = factor(N, levels = N_levels))

# Sequential ramp with the pale end trimmed, so the 10-animal series reads on white
size_pal <- setNames(colorRampPalette(brewer.pal(9, "Blues"))(7)[3:7], as.character(N_levels))
r0lab <- function(x) sprintf("R0 = %s", x)
XLIM <- c(0, 50)

# --- exported plotted data ----------------------------------------------------
panelA <- tr %>% select(R0, N, yr, remaining) %>% arrange(R0, N, yr)
panelB <- tr %>% select(R0, N, yr, p_free_vx) %>% arrange(R0, N, yr)
panelC <- sm %>% mutate(reached_half_free = !is.na(yr_half_free)) %>%
  select(R0, N, yr_half_of_fall, yr_half_free, reached_half_free, plateau) %>% arrange(R0, N)
write.csv(panelA, "output/tables/fig38_panelA_data.csv", row.names = FALSE)
write.csv(panelB, "output/tables/fig38_panelB_data.csv", row.names = FALSE)
write.csv(panelC, "output/tables/fig38_panelC_data.csv", row.names = FALSE)

thm <- function() theme_pub(11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold", size = 10),
        panel.spacing.x = unit(2.2, "lines"),
        plot.margin = margin(4, 34, 4, 16))

# --- panel A: burden remaining ------------------------------------------------
# labels go in the facet where the series separate: at R0 = 1.5 the five
# curves converge and any label set there collides whatever the repel settings
endA <- tr %>% filter(R0 == 2) %>% group_by(Nf) %>% slice_max(yr, n = 1) %>% ungroup()
annA <- data.frame(
  R0 = c(5, 3), x = c(24, 20), y = c(0.88, 0.86),
  xend = c(42, 40), yend = c(0.695, 0.60),
  lab = c("endemic plateau,\nabout 68% remains", "plateau about 58%"))
pA <- ggplot(tr, aes(yr, remaining, colour = Nf)) +
  geom_line(linewidth = 0.9) +
  geom_segment(data = annA, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey45", linewidth = 0.35) +
  geom_text(data = annA, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            colour = "grey30", size = 2.7, hjust = 0, vjust = 0, lineheight = 0.95) +
  geom_text_repel(data = endA, aes(label = N), size = 3, direction = "y", hjust = 0,
                  nudge_x = 3, segment.size = 0.3, segment.colour = "grey65",
                  min.segment.length = 0, box.padding = 0.12, xlim = c(51, NA)) +
  facet_grid(cols = vars(R0), labeller = as_labeller(r0lab)) +
  scale_colour_manual(values = size_pal) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 20, 40)) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = NULL, y = "Infection remaining") + thm()

# --- panel B: herds free ------------------------------------------------------
endB <- tr %>% filter(R0 == min(R0)) %>% group_by(Nf) %>% slice_max(yr, n = 1) %>% ungroup()
annB <- data.frame(R0 = 5, x = 3, y = 0.86, xend = 44, yend = 0.60,
                   lab = "10-animal herds stall below 60%")
pB <- ggplot(tr, aes(yr, p_free_vx, colour = Nf)) +
  geom_line(linewidth = 0.9) +
  geom_segment(data = annB, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey45", linewidth = 0.35) +
  geom_text(data = annB, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            colour = "grey30", size = 2.7, hjust = 0, vjust = -0.4) +
  geom_text_repel(data = endB, aes(label = N), size = 3, direction = "y", hjust = 0,
                  nudge_x = 3, segment.size = 0.3, segment.colour = "grey65",
                  min.segment.length = 0, box.padding = 0.12, xlim = c(51, NA)) +
  facet_grid(cols = vars(R0), labeller = as_labeller(r0lab)) +
  scale_colour_manual(values = size_pal) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 20, 40)) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = NULL, y = "Herds free of infection") +
  thm() + theme(strip.text = element_blank())

# --- panel C: milestone timings ----------------------------------------------
gc_ <- sm %>% mutate(never = is.na(yr_half_free),
                     xend = ifelse(never, 50, yr_half_free))
decade <- median(sm$yr_half_of_fall, na.rm = TRUE)
pC <- ggplot(gc_, aes(y = Nf, colour = Nf)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.16) +
  geom_vline(xintercept = decade, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_segment(aes(x = yr_half_of_fall, xend = xend, yend = Nf), linewidth = 1.1) +
  geom_segment(data = filter(gc_, never), aes(x = 50, xend = 52.2, yend = Nf),
               arrow = arrow(length = unit(0.09, "cm"), type = "closed"), linewidth = 1.1) +
  geom_text(data = filter(gc_, never), aes(x = 53.2, y = Nf, label = "> 50 y"),
            hjust = 0, size = 2.5, colour = "grey35") +
  geom_point(aes(x = yr_half_of_fall), size = 2.4) +
  geom_point(data = filter(gc_, !never), aes(x = xend), shape = 21, fill = "white",
             size = 2.4, stroke = 1.05) +
  geom_text(data = data.frame(R0 = min(sm$R0)), aes(x = 0.8, y = 5.40, label = "first decade"),
            inherit.aes = FALSE, hjust = 0, size = 2.6, colour = "grey40") +
  facet_grid(cols = vars(R0), labeller = as_labeller(r0lab)) +
  scale_colour_manual(values = size_pal) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50)) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = "Years since the start of vaccination", y = "Herd size (animals)") +
  thm() + theme(strip.text = element_blank(), panel.grid.major.y = element_blank(),
                panel.spacing.x = unit(2.9, "lines"), plot.margin = margin(4, 46, 4, 16))

# --- in-plot key for panel C --------------------------------------------------
key <- ggplot() +
  geom_point(aes(x = 0.5, y = 1), size = 2.4, colour = "#105BA4") +
  annotate("text", x = 1.2, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "half the burden reduction achieved") +
  geom_point(aes(x = 12.5, y = 1), shape = 21, fill = "white", size = 2.4, stroke = 1.05, colour = "#105BA4") +
  annotate("text", x = 13.2, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "half of herds free of infection") +
  annotate("segment", x = 23.4, xend = 24.6, y = 1, yend = 1, colour = "#105BA4",
           linewidth = 1.1, arrow = arrow(length = unit(0.09, "cm"), type = "closed")) +
  annotate("text", x = 25.3, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "not reached within 50 years") +
  scale_x_continuous(limits = c(0, 40)) + scale_y_continuous(limits = c(0.9, 1.1)) +
  theme_void() + theme(plot.margin = margin(0, 34, 2, 16))

title <- ggdraw() + draw_label("Benefit is banked within a decade in every stratum",
                               fontface = "bold", size = 13, x = 0.012, hjust = 0, vjust = 0.5)
body <- plot_grid(pA, pB, pC, ncol = 1, labels = c("A", "B", "C"), label_size = 14,
                  label_x = 0.002, rel_heights = c(1, 1, 1.05), align = "v", axis = "lr")
fig <- plot_grid(title, body, key, ncol = 1, rel_heights = c(0.05, 1, 0.045))
ggsave("output/figs/fig38_strata_detail.png", fig, width = 11, height = 11, dpi = 200, bg = "white")
cat("fig38 written; panel data exported to output/tables/fig38_panel[ABC]_data.csv\n")
