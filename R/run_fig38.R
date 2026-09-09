# Figure 38 (revised): the twenty demographic strata of the benefit-timing
# analysis. Plotting layer only - reads the tables written by
# R/run_benefit_timing.R. Run from the project root: Rscript R/run_fig38.R
#
# Revisions in this version:
#   - title claim matched to the data (halfway point spans 2.5-11.5 years)
#   - herd-size ramp on an even OKLCH lightness scale (greyscale-safe)
#   - R0 across columns in all three panels; strip labels parsed as R[0]
#   - series direct-labelled in the R0 = 2 facet only, dark text + coloured leader
#   - panel C: first-decade band, no competing median line, one ink, flat
#     segment ends with "> 50 y", concentric marker for same-year strata
#   - unified 0-50 x breaks; A's redundant x tick labels dropped

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(cowplot); library(ggrepel)
})
source("R/root.R")
source("R/model.R")   # theme_pub()

N_levels <- c(10, 20, 50, 100, 200)
tr <- read.csv("output/tables/benefit_timing.csv")         %>% mutate(Nf = factor(N, levels = N_levels))
sm <- read.csv("output/tables/benefit_timing_summary.csv") %>% mutate(Nf = factor(N, levels = N_levels))

# Herd-size ramp: one steel hue on an even OKLCH lightness scale (L 0.88 -> 0.38).
# Monotone lightness means the order survives greyscale printing.
size_pal <- setNames(c("#BADCFA", "#7CB6E4", "#4E8DCC", "#3864A5", "#2A4072"),
                     as.character(N_levels))
# Equivalent generator, if you would rather derive than hard-code:
#   size_pal <- setNames(
#     farver::encode_colour(farver::convert_colour(
#       cbind(l = c(0.88, 0.755, 0.63, 0.505, 0.38),
#             c = c(0.055, 0.090, 0.115, 0.115, 0.090),
#             h = c(245,   243,   250,   258,   265)),
#       from = "oklch", to = "rgb"), from = "rgb"),
#     as.character(N_levels))

INK   <- "#3864A5"   # single ink for panel C: colour there is redundant with the
                     # labelled rows, and the pale steps fall below 3:1 at mark size
XLIM  <- c(0, 50)
BRK   <- c(0, 10, 20, 30, 40, 50)
LABF  <- 2           # facet in which the series are direct-labelled

thm <- function() theme_pub(11) +
  theme(legend.position   = "none",
        panel.grid.minor  = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_text(hjust = 0, face = "bold", size = 10),
        panel.spacing.x   = unit(2.4, "lines"),
        plot.margin       = margin(4, 52, 4, 16))

facet_r0 <- facet_grid(cols = vars(R0), labeller = label_bquote(cols = R[0] == .(R0)))

# ---- panels A and B ---------------------------------------------------------
ends <- function(d, yvar) d %>% filter(R0 == LABF) %>% group_by(Nf) %>%
  slice_max(yr, n = 1) %>% ungroup() %>% rename(yval = {{ yvar }})

lines_panel <- function(yvar, ylab, ann, xaxis) {
  e <- ends(tr, {{ yvar }})
  ggplot(tr, aes(yr, {{ yvar }}, colour = Nf)) +
    geom_line(linewidth = 1) +
    geom_segment(data = ann, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "grey45", linewidth = 0.35) +
    geom_text(data = ann, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
              colour = "grey30", size = 2.7, hjust = 0, vjust = 0, lineheight = 0.95) +
    # dark label text, coloured leader: the leader carries the curve association
    geom_text_repel(data = e, aes(label = N, segment.colour = Nf), y = e$yval,
                    colour = "grey10", size = 3, fontface = "bold",
                    direction = "y", hjust = 0, nudge_x = 4, xlim = c(51, NA),
                    ylim = c(0, 1), force = 2, max.overlaps = Inf,
                    box.padding = 0.12, segment.size = 0.35, min.segment.length = 0) +
    facet_r0 +
    scale_colour_manual(values = size_pal, aesthetics = c("colour", "segment.colour")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_x_continuous(breaks = BRK) +
    coord_cartesian(xlim = XLIM, clip = "off") +
    labs(x = NULL, y = ylab) + thm() +
    (if (xaxis) theme() else theme(axis.text.x = element_blank()))
}

annA <- data.frame(R0 = c(3, 5), x = c(20, 18), y = c(0.86, 0.88),
                   xend = c(40, 42), yend = c(0.576, 0.674),
                   lab = c("herds \u2265 50 animals\nplateau at 53-58%",
                           "\u2265 50 animals: endemic\nplateau, 66-67% remains"))
annB <- data.frame(R0 = 5, x = 26, y = 0.80, xend = 45, yend = 0.53,
                   lab = "10-animal herds\nstall below 60%")

pA <- lines_panel(remaining,  "Infection remaining",     annA, xaxis = FALSE)
pB <- lines_panel(p_free_vx,  "Herds free of infection", annB, xaxis = TRUE) +
  theme(strip.text = element_blank())

# ---- panel C: milestone timings ---------------------------------------------
# some strata reach both milestones in the same year (and the second can land
# marginally earlier): one concentric marker, never a backwards segment
cc <- sm %>% mutate(
  never = is.na(yr_half_free),
  same  = !never & abs(yr_half_free - yr_half_of_fall) < 0.5,
  xend  = ifelse(never, 50, pmax(yr_half_free, yr_half_of_fall)))

pC <- ggplot(cc, aes(y = Nf)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.16) +
  geom_segment(aes(x = yr_half_of_fall, xend = xend, yend = Nf),
               colour = INK, linewidth = 1.1) +
  geom_point(data = filter(cc, !never, !same), aes(x = xend),
             shape = 21, fill = "white", colour = INK, size = 2.4, stroke = 1.05) +
  geom_point(data = filter(cc, !same), aes(x = yr_half_of_fall),
             colour = INK, size = 2.4) +
  geom_point(data = filter(cc, same), aes(x = yr_half_of_fall), colour = INK, size = 3.3) +
  geom_point(data = filter(cc, same), aes(x = yr_half_of_fall), colour = "white", size = 1.2) +
  geom_text(data = filter(cc, never), aes(x = 51.5, label = "> 50 y"),
            hjust = 0, size = 2.5, colour = "grey35") +
  geom_text(data = data.frame(R0 = 1.5), aes(x = 0.6, y = 5.42, label = "first decade"),
            inherit.aes = FALSE, hjust = 0, size = 2.6, colour = "grey40") +
  facet_r0 +
  scale_x_continuous(breaks = BRK) +
  coord_cartesian(xlim = XLIM, clip = "off") +
  labs(x = "Years since the start of vaccination", y = "Herd size (animals)") +
  thm() + theme(panel.grid.major.y = element_blank(),
                plot.margin = margin(4, 58, 4, 16))

# ---- key, title, assembly ---------------------------------------------------
# One key row, for panel C's markers only. The herd-size ramp swatch is gone:
# panels A and B are directly labelled at R0 = 2, so a colour legend repeats
# information the reader already has. The former footer is gone too, since it
# restated these same four marks in prose.
key <- ggplot() +
  geom_point(aes(0.4, 1), size = 2.4, colour = INK) +
  annotate("text", x = 1.2, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "half the burden reduction achieved") +
  geom_point(aes(13.0, 1), shape = 21, fill = "white", size = 2.4, stroke = 1.05, colour = INK) +
  annotate("text", x = 13.8, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "half of herds free of infection") +
  geom_point(aes(24.5, 1), size = 3.3, colour = INK) +
  geom_point(aes(24.5, 1), size = 1.2, colour = "white") +
  annotate("text", x = 25.3, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "both milestones in the same year") +
  annotate("text", x = 36.5, y = 1, hjust = 0, size = 2.6, colour = "grey40", label = "> 50 y") +
  annotate("text", x = 39.1, y = 1, hjust = 0, size = 2.9, colour = "grey25",
           label = "not reached within 50 years") +
  scale_x_continuous(limits = c(0, 60)) + scale_y_continuous(limits = c(0.9, 1.1)) +
  theme_void() + theme(legend.position = "none", plot.margin = margin(0, 16, 0, 16))

head_ <- ggdraw() +
  draw_label("Half the burden reduction arrives within roughly a decade in every stratum",
             fontface = "bold", size = 13.5, x = 0.012, y = 0.66, hjust = 0) +
  draw_label(paste("Across the twenty strata the halfway point falls between 2.5 and 11.5 years;",
                   "the shaded band in C is that first decade."),
             size = 9, colour = "grey30", x = 0.012, y = 0.22, hjust = 0)

# named panels, not body: base R already defines body() and shadowing it here
# silently yields a function object and an empty figure
panels <- plot_grid(pA, pB, pC, ncol = 1, labels = c("A", "B", "C"),
                    label_size = 14, label_x = 0.002,
                    rel_heights = c(1, 1, 1.05), align = "v", axis = "lr")

fig <- plot_grid(head_, panels, key, ncol = 1,
                 rel_heights = c(0.065, 1, 0.042))

# Plotted data for each panel, as required alongside the figure
write.csv(tr %>% select(R0, N, yr, remaining) %>% arrange(R0, N, yr),
          "output/tables/fig38_panelA_data.csv", row.names = FALSE)
write.csv(tr %>% select(R0, N, yr, p_free_vx) %>% arrange(R0, N, yr),
          "output/tables/fig38_panelB_data.csv", row.names = FALSE)
write.csv(sm %>% mutate(reached_half_free = !is.na(yr_half_free)) %>%
            select(R0, N, yr_half_of_fall, yr_half_free, reached_half_free, plateau) %>%
            arrange(R0, N),
          "output/tables/fig38_panelC_data.csv", row.names = FALSE)

ggsave("output/figs/fig38_strata_detail.png", fig,
       width = 11.8, height = 10.2, dpi = 300, bg = "white")
ggsave("output/figs/fig38_strata_detail.pdf", fig,
       width = 11.8, height = 10.2, device = cairo_pdf, bg = "white")
cat("fig38 written (png + pdf)\n")

# Acceptance check: convert the PNG to greyscale and confirm all five herd sizes
# stay distinguishable in the R0 = 1.5 facet of panel A - that convergent bundle
# is the failure mode the new ramp fixes.
#   magick convert output/figs/fig38_strata_detail.png -colorspace Gray /tmp/fig38_grey.png

# Notes
# - aes(segment.colour = Nf) needs ggrepel >= 0.9.0; on older versions drop it
#   from aes() and pass segment.colour = "grey60" as a plain argument.
# - theme_pub() comes from R/root.R; thm() assumes it does not set plot.margin.
