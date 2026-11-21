# =============================================================================
# Sticker motifs for Tweedieverse.
# One model, many distributions: as the Tweedie power p varies, the family
# covers Normal, Poisson, compound Poisson-gamma, Gamma and inverse Gaussian.
# =============================================================================
library(ggplot2)

TW_GOLD <- "#F4B942"

# Simulated Tweedie: P(Y = 0) and a smooth density for Y > 0 (scaled to P(Y > 0))
.tweedie_parts <- function(p, mu = 1, phi = 1, n = 4e5, seed = 1, xmax = 4) {
  set.seed(seed)
  lambda <- mu^(2 - p) / (phi * (2 - p))
  shape  <- (2 - p) / (p - 1)
  scale  <- phi * (p - 1) * mu^(p - 1)
  N <- rpois(n, lambda)
  y <- numeric(n); pos <- N > 0
  y[pos] <- rgamma(sum(pos), shape = N[pos] * shape, scale = scale)
  d <- density(y[pos], from = 0, to = xmax, n = 300, bw = 0.12)
  list(p0 = exp(-lambda), dens = data.frame(x = d$x, y = d$y * mean(pos)))
}

.finish_tw <- function(p, xlim, ylim) {
  p + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void() + hexSticker::theme_transparent()
}

# -----------------------------------------------------------------------------
# One family, many distributions: the Tweedie power p spans
#   p = 0 Normal | p = 1 Poisson | 1 < p < 2 compound Poisson-gamma |
#   p = 2 Gamma  | p = 3 inverse Gaussian
# Each member is drawn as a small "icon" of its distribution shape.
# -----------------------------------------------------------------------------
TW_MEMBERS <- data.frame(
  name = c("normal", "poisson", "cpg", "gamma", "invgauss"),
  p    = c(0, 1, 1.5, 2, 3),
  col  = c("#4EA8DE", "#2EC4B6", "#F4B942", "#FF7A7A", "#9B8CFF"),
  stringsAsFactors = FALSE
)

# layers drawing one distribution icon inside the box centred at (cx, cy)
.dist_icon <- function(name, cx, cy, w = 1.2, h = 0.8, col, lw = 1.1) {
  base <- data.frame(x = cx - w / 2, xend = cx + w / 2, y = cy - h / 2)
  X <- function(u) cx - w / 2 + u * w           # u in [0, 1] -> x
  Y <- function(v) cy - h / 2 + v * h           # v in [0, 1] -> y
  u <- seq(0, 1, length.out = 120)
  curve_layer <- function(v, fill = TRUE) {
    d <- data.frame(x = X(u), y = Y(v / max(v) * 0.92))
    d$y0 <- Y(0)
    list(if (fill) geom_ribbon(data = d, aes(x, ymin = y0, ymax = y), fill = col, alpha = 0.22),
         geom_line(data = d, aes(x, y), colour = col, linewidth = lw, lineend = "round"))
  }
  out <- switch(name,
    normal   = curve_layer(dnorm(u, 0.5, 0.15)),
    gamma    = curve_layer(dgamma(u * 6, shape = 2.2, rate = 1.3)),
    invgauss = {
      t <- pmax(u * 4, 1e-3); m <- 1; lam <- 1.2
      curve_layer(sqrt(lam / (2 * pi * t^3)) * exp(-lam * (t - m)^2 / (2 * m^2 * t)))
    },
    poisson  = {
      k <- 0:7; v <- dpois(k, 2.6); v <- v / max(v) * 0.92
      bw <- w / 8 * 0.62
      d <- data.frame(xmin = X((k + 0.5) / 8) - bw / 2, xmax = X((k + 0.5) / 8) + bw / 2,
                      ymin = Y(0), ymax = Y(v))
      list(geom_rect(data = d, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     fill = col))
    },
    cpg = {
      tw <- .tweedie_parts(1.5, xmax = 3.2)
      v <- tw$dens$y; v <- v / max(v) * 0.7
      d <- data.frame(x = X(0.14 + 0.86 * tw$dens$x / 3.2), y = Y(v), y0 = Y(0))
      sp <- data.frame(x = X(0.04), y = Y(0.92))
      list(geom_ribbon(data = d, aes(x, ymin = y0, ymax = y), fill = col, alpha = 0.22),
           geom_line(data = d, aes(x, y), colour = col, linewidth = lw, lineend = "round"),
           geom_segment(data = sp, aes(x = x, xend = x, y = Y(0), yend = y),
                        colour = col, linewidth = lw * 1.1),
           geom_point(data = sp, aes(x, y), colour = col, size = lw * 1.5))
    })
  c(list(geom_segment(data = base, aes(x = x, xend = xend, y = y, yend = y),
                      colour = "white", alpha = 0.55, linewidth = 0.35)), out)
}

# A. Spectrum: the five family members on an arc above one p-axis; each sits
#    over its own value of the Tweedie power p (dots placed at their true p)
motif_family_arc <- function(m = TW_MEMBERS, family = "sticker_font",
                             labels = c("0", "1", "1\u20132", "2", "3"), label_size = 4.8) {
  n <- nrow(m); ang <- seq(165, 15, length.out = n) * pi / 180
  R <- 2.45; cx <- cos(ang) * R; cy <- sin(ang) * R * 0.72 + 0.55
  p <- ggplot()
  for (i in seq_len(n)) p <- p + .dist_icon(m$name[i], cx[i], cy[i], w = 1.0, h = 0.66,
                                            col = m$col[i])
  ax_y <- -0.6
  ax <- data.frame(x = -2.2 + m$p / 3 * 4.4, y = ax_y, lab = labels, col = m$col)
  links <- data.frame(x = ax$x, y = ax$y, xend = cx, yend = cy - 0.38)
  p <- p +
    geom_segment(aes(x = -2.45, xend = 2.45, y = ax_y, yend = ax_y),
                 colour = "white", linewidth = 0.8) +
    geom_segment(data = links, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "white", alpha = 0.3, linewidth = 0.4, linetype = "22") +
    geom_point(data = ax, aes(x, y), colour = ax$col, size = 2.6) +
    geom_text(data = ax, aes(x, y - 0.42, label = lab), colour = "white",
              family = family, fontface = "bold", size = label_size) +
    annotate("text", x = 2.78, y = ax_y, label = "p", colour = "white",
             family = family, fontface = "bold.italic", size = label_size * 1.15)
  .finish_tw(p, c(-3.0, 3.1), c(-1.3, max(cy) + 0.42))
}
