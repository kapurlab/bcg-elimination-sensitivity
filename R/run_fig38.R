# Figure 38: the 20 demographic strata of the benefit-timing analysis, drawn
# separately by herd size and within-herd R0.
#
# Plotting layer only. The model output is read from the tables written by
# R/run_benefit_timing.R and is not recomputed here, so this script can be
# iterated on freely without touching the simulation.
#
# Run from the project root: Rscript R/run_fig38.R (seconds)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(cowplot); library(ggrepel)
})
source("R/root.R"); source("R/model.R")

N_levels <- c(10, 20, 50, 100, 200)
tr <- read.csv("output/tables/benefit_timing.csv") %>% mutate(Nf = factor(N, levels = N_levels))
sm <- read.csv("output/tables/benefit_timing_summary.csv") %>% mutate(Nf = factor(N, levels = N_levels))

# Herd-size ramp: single steel hue, five steps on an even OKLCH lightness
# scale (L = 0.88 -> 0.38). Monotone lightness means the order survives
# greyscale printing; wider range than the trimmed Brewer subset it replaces.
size_pal <- setNames(
  c("#BADCFA", "#7CB6E4", "#4E8DCC", "#3864A5", "#2A4072"),
  as.character(N_levels))   # 10, 20, 50, 100, 200

r0_strip <- as_labeller(function(x) sprintf("R[0] == %s", x), default = label_parsed)
XLIM <- c(0, 50)
XBREAKS <- c(0, 10, 20, 30, 40, 50)

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
        strip.text = element_text(hjust = 0, size = 10),
        panel.spacing.x = unit(3.0, "lines"),
        plot.margin = margin(4, 52, 4, 16))

# Direct labels sit in the R0 = 2 facet of both A and B: the five series
# separate there in each panel, whereas at R0 = 1.5 the burden curves converge
# to within a few points and any label set collides whatever the repel settings.
LABEL_FACET <- 2
ends <- tr %>% filter(R0 == LABEL_FACET) %>% group_by(Nf) %>% slice_max(yr, n = 1) %>% ungroup()
repel_layer <- function() {
  geom_text_repel(data = ends, aes(label = N), size = 3, direction = "y", hjust = 0,
                  nudge_x = 5, force = 2, max.overlaps = Inf, ylim = c(0, 1),
                  xlim = c(51, NA), segment.size = 0.3, segment.colour = "grey65",
                  min.segment.length = 0, box.padding = 0.12)
}

# --- panel A: burden remaining ------------------------------------------------
# Plateau values are quoted as the range across herds of 50 or more, which is
# what the flat bundle in each facet actually spans (53 to 58% at R0 = 3,
# 66 to 67% at R0 = 5). Labelling each curve at its terminus was considered and
# rejected: the three largest herds converge to within one point at R0 = 5, so
# per-curve labels would overlap exactly where the plateau matters.
annA <- data.frame(
  R0 = c(5, 3), x = c(19, 19), y = c(0.88, 0.84),
  xend = c(42, 40), yend = c(0.70, 0.60),
  lab = c("50+ animals: endemic\nplateau at 66 to 67%", "50+ animals plateau\nat 53 to 58%"))
pA <- ggplot(tr, aes(yr, remaining, colour = Nf)) +
  geom_line(linewidth = 1.0) +
  geom_segment(data = annA, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey45", linewidth = 0.35) +
  geom_text(data = annA, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            colour = "grey30", size = 2.7, hjust = 0, vjust = 0, lineheight = 0.95) +
  repel_layer() +
  facet_grid(cols = vars(R0), labeller = r0_strip) +
  scale_colour_manual(values = size_pal) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = XBREAKS) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = NULL, y = "Infection remaining") +
  thm() + theme(axis.text.x = element_blank())

# --- panel B: herds free ------------------------------------------------------
annB <- data.frame(R0 = 5, x = 24, y = 0.76, xend = 47, yend = 0.61,
                   lab = "10-animal herds\nstall below 60%")
pB <- ggplot(tr, aes(yr, p_free_vx, colour = Nf)) +
  geom_line(linewidth = 1.0) +
  geom_segment(data = annB, aes(x = x, y = y, xend = xend, yend = yend),
               inherit.aes = FALSE, colour = "grey45", linewidth = 0.35) +
  geom_text(data = annB, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            colour = "grey30", size = 2.7, hjust = 0, vjust = 0, lineheight = 0.95) +
  repel_layer() +
  facet_grid(cols = vars(R0), labeller = r0_strip) +
  scale_colour_manual(values = size_pal) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(breaks = XBREAKS) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = NULL, y = "Herds free of infection") +
  thm() + theme(strip.text = element_blank())

# --- panel C: milestone timings ----------------------------------------------
# "Not reached" is encoded once, by the flat segment end at the axis limit plus
# the "> 50 y" label. No arrowhead, and no missing symbol to notice.
gc_ <- sm %>% mutate(never = is.na(yr_half_free), xend = ifelse(never, 50, yr_half_free))
pC <- ggplot(gc_, aes(y = Nf, colour = Nf)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.16) +
  geom_segment(aes(x = yr_half_of_fall, xend = xend, yend = Nf), linewidth = 1.1) +
  geom_text(data = filter(gc_, never), aes(x = 51.5, y = Nf, label = "> 50 y"),
            hjust = 0, size = 2.5, colour = "grey35") +
  geom_point(aes(x = yr_half_of_fall), size = 2.4) +
  geom_point(data = filter(gc_, !never), aes(x = xend), shape = 21, fill = "white",
             size = 2.4, stroke = 1.05) +
  geom_text(data = data.frame(R0 = min(sm$R0)), aes(x = 9.4, y = 5.35, label = "first decade"),
            inherit.aes = FALSE, hjust = 1, size = 2.6, colour = "grey40") +
  facet_grid(cols = vars(R0), labeller = r0_strip) +
  scale_colour_manual(values = size_pal) +
  scale_x_continuous(breaks = XBREAKS) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = "Years since the start of vaccination", y = "Herd size (animals)") +
  thm() + theme(panel.grid.major.y = element_blank())

# --- in-plot key for panel C --------------------------------------------------
key_col <- unname(size_pal["100"])
key <- ggplot() +
  geom_point(aes(x = 0.5, y = 1), size = 2.4, colour = key_col) +
  annotate("text", x = 1.3, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "half the burden reduction achieved") +
  geom_point(aes(x = 13.5, y = 1), shape = 21, fill = "white", size = 2.4,
             stroke = 1.05, colour = key_col) +
  annotate("text", x = 14.3, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "half of herds free of infection") +
  annotate("text", x = 24.5, y = 1, hjust = 0, size = 2.9, colour = "grey35", label = "> 50 y") +
  annotate("text", x = 25.85, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = ", not reached within 50 years") +
  scale_x_continuous(limits = c(0, 40)) + scale_y_continuous(limits = c(0.9, 1.1)) +
  theme_void() + theme(plot.margin = margin(0, 52, 2, 16))

head_row <- ggdraw() +
  draw_label("Half the burden reduction arrives within roughly a decade in every stratum",
             fontface = "bold", size = 13, x = 0.012, hjust = 0, y = 0.70) +
  draw_label(paste("Across the 20 strata that milestone falls between 2.5 and 11.5 years.",
                   "Freedom from infection comes far later, or not at all."),
             size = 9.5, colour = "grey30", x = 0.012, hjust = 0, y = 0.26)
body <- plot_grid(pA, pB, pC, ncol = 1, labels = c("A", "B", "C"), label_size = 14,
                  label_x = 0.002, rel_heights = c(1, 1, 1.05), align = "v", axis = "lr")
fig <- plot_grid(head_row, body, key, ncol = 1, rel_heights = c(0.075, 1, 0.045))
ggsave("output/figs/fig38_strata_detail.png", fig, width = 11, height = 11, dpi = 200, bg = "white")
cat("fig38 written; panel data exported to output/tables/fig38_panel[ABC]_data.csv\n")
