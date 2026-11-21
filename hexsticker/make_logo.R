# =============================================================================
# Tweedieverse hex logo: one Tweedie model, many distributions along the power p
#   Run from the package root (RStudio project):  source("hexsticker/make_logo.R")
# Requirements: install.packages(c("hexSticker", "showtext", "ggplot2", "usethis", "pkgdown"))
# =============================================================================
source("hexsticker/make_hexsticker.R")
source("hexsticker/motifs_tw.R")

make_pkg_sticker(
  package  = "Tweedieverse",
  subplot  = motif_family_arc(),     # Normal, Poisson, CPG, Gamma, inverse Gaussian
  s_width  = 1.4, s_height = 0.95, s_y = 0.82,
  filename = "hexsticker/logo.png"
)

# install it as the package / website logo and refresh the browser-tab icons
usethis::use_logo("hexsticker/logo.png", retina = TRUE)
pkgdown::build_favicons(overwrite = TRUE)
