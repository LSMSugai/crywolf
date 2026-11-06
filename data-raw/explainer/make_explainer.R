# Figures of vignettes/articles/how-crywolf-works.Rmd, in the palette and
# theme of the plots of the package. Needs ggplot2, patchwork and Graphviz
# (`dot`). Run from the package root:
#   Rscript data-raw/explainer/make_explainer.R
library(ggplot2)
library(patchwork)

ink <- "#241042"
pink <- "#FF3EA5"
cyan <- "#2DE2E6"
none <- "grey91"
out <- file.path("vignettes", "articles", "figures")

# the values of the worked examples
p_occupied <- function(n1, n, psi = 0.5, p11 = 0.6, p10 = 0.05) {
  occupied <- psi * p11^n1 * (1 - p11)^(n - n1)
  empty <- (1 - psi) * p10^n1 * (1 - p10)^(n - n1)
  occupied / (occupied + empty)
}

base <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title.position = "plot",
        plot.title = element_text(size = 12, face = "bold", colour = ink))
bare <- base + theme(panel.grid = element_blank(), axis.text = element_blank(),
                     axis.title = element_blank())
save <- function(p, name, width, height) {
  ggsave(file.path(out, paste0("how-", name, ".png")), p, width = width,
         height = height, dpi = 200, bg = "white")
}
wide <- 8.2
narrow <- 5

# the problem: six sites over ten nights
hits <- list(S1 = c(1, 3, 4, 6, 8, 9), S2 = c(2, 5, 8), S3 = integer(0),
             S4 = c(3, 8), S5 = 7, S6 = integer(0))
occupied <- c(S1 = TRUE, S2 = TRUE, S3 = TRUE, S4 = FALSE, S5 = FALSE, S6 = FALSE)
nights <- expand.grid(night = 1:10, site = names(hits), stringsAsFactors = FALSE)
nights$state <- mapply(function(n, s) {
  if (!n %in% hits[[s]]) "no detection"
  else if (occupied[[s]]) "true detection"
  else "false positive"
}, nights$night, nights$site)
nights$state <- factor(nights$state, c("true detection", "false positive", "no detection"))
sites <- data.frame(site = names(hits), n1 = lengths(hits),
                    truth = ifelse(occupied, "occupied", "empty"))
sites$naive <- ifelse(sites$n1 > 0, "occupied", "empty")
sites$p <- sprintf("%.2f", p_occupied(sites$n1, 10))
p <- ggplot(nights, aes(night, site, fill = state)) +
  geom_tile(width = 0.9, height = 0.78) +
  geom_text(data = sites, aes(-1.1, site, label = truth), inherit.aes = FALSE,
            size = 3.6, colour = ink) +
  geom_text(data = sites, aes(12.2, site, label = naive), inherit.aes = FALSE,
            size = 3.6, colour = ink) +
  geom_text(data = sites, aes(14.6, site, label = p), inherit.aes = FALSE,
            size = 3.6, fontface = "bold", colour = ink) +
  annotate("text", x = c(-1.1, 5.5, 12.2, 14.6), y = 6.75, size = 3.4,
           colour = "grey30", lineheight = 0.9, vjust = 0,
           label = c("truth", "ten nights", "one detection\nis enough", "P(occupied)")) +
  scale_y_discrete(limits = rev(names(hits))) +
  scale_x_continuous(limits = c(-2.2, 15.6), expand = c(0, 0)) +
  scale_fill_manual(values = c("true detection" = ink, "false positive" = pink,
                               "no detection" = none), name = NULL) +
  coord_cartesian(ylim = c(0.5, 7.5), clip = "off") +
  bare + theme(axis.text.y = element_text(), legend.position = "bottom")
save(p, "problem", wide, 3.6)

# the mixture: expected sites out of 200. p11 is lower here than in the
# worked examples, so that the two groups overlap
k <- 0:15
groups <- rbind(data.frame(k, n = 100 * dbinom(k, 20, 0.05), group = "empty sites"),
                data.frame(k, n = 100 * dbinom(k, 20, 0.25), group = "occupied sites"))
p <- ggplot(groups, aes(k, n, fill = group)) +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = cyan, alpha = 0.18) +
  geom_col(position = position_dodge(width = 0.86), width = 0.8) +
  annotate("text", x = 2.5, y = 42, label = "could be either", size = 3.6, colour = ink) +
  scale_fill_manual(values = c("empty sites" = pink, "occupied sites" = ink), name = NULL) +
  scale_x_continuous(breaks = seq(0, 15, 5)) +
  labs(x = "nights with a detection, out of 20", y = "expected sites, of 200") +
  base + theme(legend.position = "top", legend.justification = "left",
               panel.grid.major.x = element_blank())
save(p, "mixture", narrow, 3.2)

# P(occupied) against the number of detections
curve <- data.frame(n1 = 0:10)
curve$p <- p_occupied(curve$n1, 10)
worked <- curve[curve$n1 %in% 2:3, ]
worked$label <- sprintf("%d of 10: %.2f", worked$n1, worked$p)
p <- ggplot(curve, aes(n1, p)) +
  geom_line(colour = ink, linewidth = 0.9) +
  geom_point(colour = ink, size = 2) +
  geom_point(data = worked, colour = pink, size = 3.4) +
  geom_text(data = worked, aes(label = label), hjust = 0, nudge_x = 0.35,
            nudge_y = c(0, -0.06), size = 3.6, colour = ink) +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "nights with a detection, out of 10", y = "P(occupied)") +
  base + theme(panel.grid.major.x = element_blank())
save(p, "p-occupied", narrow, 3)

# the two designs: verified occasions are occasions of their own
fills <- c("0" = none, "1" = ink, "F" = pink, "v0" = cyan, "v1" = ink)
multiple <- rbind(
  data.frame(x = 1:10, row = "automated",
             f = c("0", "1", "0", "0", "F", "0", "0", "1", "0", "0")),
  data.frame(x = c(4, 6, 8), row = "verified", f = c("v0", "v1", "v0")))
multiple$row <- factor(multiple$row, c("verified", "automated"))
top <- ggplot(multiple, aes(x, row, fill = f)) +
  geom_tile(width = 0.9, height = 0.7) +
  annotate("text", x = 9.2, y = 1, hjust = 0, size = 3.5, colour = ink, lineheight = 0.95,
           label = "dark: heard, so z = 1\ncyan: nothing heard") +
  scale_fill_manual(values = fills, guide = "none") +
  scale_x_continuous(limits = c(0.4, 14)) +
  labs(title = "Multiple methods") +
  bare + theme(axis.text.y = element_text())
states <- data.frame(x = 0:3 * 3.2, f = c("0", "F", "1", "1"))
bottom <- ggplot() +
  geom_tile(data = states, aes(x, 1, fill = f), width = 1.6, height = 0.8) +
  annotate("rect", xmin = 9.6, xmax = 10.4, ymin = 0.6, ymax = 1.4, fill = pink) +
  annotate("text", x = states$x, y = 0.35, size = 3.5, colour = ink, lineheight = 0.9,
           label = c("0\nnothing flagged", "1\nall false", "2\nall true", "3\ntrue and false")) +
  scale_fill_manual(values = fills, guide = "none") +
  scale_y_continuous(limits = c(-0.1, 1.5)) +
  labs(title = "Observation confirmation") +
  bare
save(top / bottom + plot_layout(heights = c(1, 1.1)), "designs", wide, 3.6)

# the threshold: an occasion is a detection when any prediction reaches it
scores <- data.frame(
  occasion = rep(1:5, c(3, 2, 1, 4, 0)),
  score = c(0.82, 0.52, 0.36, 0.40, 0.22, 0.68, 0.58, 0.45, 0.30, 0.12))
scores$above <- scores$score >= 0.5
y <- data.frame(occasion = 1:5, label = paste("y =", c(1, 0, 1, 1, 0)))
p <- ggplot(scores, aes(occasion, score)) +
  geom_hline(yintercept = 0.5, colour = pink, linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(colour = above), size = 2.8) +
  geom_text(data = y, aes(occasion, -0.1, label = label), size = 3.6, colour = ink) +
  annotate("text", x = 5.45, y = 0.54, label = "threshold", hjust = 1, size = 3.5, colour = ink) +
  scale_colour_manual(values = c("TRUE" = ink, "FALSE" = "grey60"), guide = "none") +
  scale_x_continuous(breaks = 1:5, limits = c(0.5, 5.5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(-0.15, 1)) +
  labs(x = "occasion", y = "score") +
  base + theme(panel.grid.major.x = element_blank())
save(p, "threshold", narrow, 3)

# the covariates: illustrative curves
slope <- data.frame(x = seq(-2, 2, length.out = 100))
slope$true <- plogis(0.2 - 1.8 * slope$x)
slope$constant <- plogis(0.2 - 0.6 * slope$x)
p <- ggplot(slope, aes(x)) +
  geom_line(aes(y = true), colour = ink, linewidth = 0.9) +
  geom_line(aes(y = constant), colour = pink, linewidth = 0.9, linetype = "dashed") +
  annotate("text", x = c(-1.05, 0.25), y = c(0.97, 0.6), hjust = 0, size = 3.6, colour = ink,
           label = c("true occupancy", "fit with constant p10")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "covariate tied to occupancy and to false alarms", y = expression(psi)) +
  base + theme(panel.grid.major.x = element_blank())
save(p, "covariates", narrow, 3)

# the diagrams: the two stories of a site, and the routes to results
for (name in c("model", "results")) {
  system2("dot", c("-Tpng", "-Gdpi=200", file.path("data-raw", "explainer", paste0(name, ".dot")),
                   "-o", file.path(out, paste0("how-", name, ".png"))))
}
