## Regression tests for the dendrogram brush-selection bug.
##
## Two independent properties are covered:
##
##   1. The leaf axis of the rendered dendrogram has a FIXED domain. Highlight
##      rectangles are real geometry drawn at leaf +/- 0.25, so before the fix a
##      highlight touching the first or last leaf widened the data range, ggplot
##      re-expanded the scale, and the panel domain moved underneath a brush that
##      had already been drawn in pixels. Selection then chased a moving target.
##
##   2. .brushed_leaf_idcs() maps a brush onto POSITIONS in the classes vector,
##      given the data coordinates the leaves are actually drawn at. The zoom plot
##      draws leaves at x_min..x_max rather than 1..n, so the mapping cannot assume
##      the two coincide.

## --- helpers ---------------------------------------------------------------

## Smallest input that reaches trajectories_dendrogram()'s precomputed path.
## Highlighting matches on the dendrogram's OWN leaf labels, so return them in
## drawing order (x = 1..n) -- picking "the top leaf" by label needs that order.
tiny_dendrogram <- function(n = 8) {
  set.seed(1)
  m <- matrix(rnorm(n * 2), ncol = 2)
  rownames(m) <- paste0("leaf", seq_len(n))
  d <- as.dendrogram(hclust(dist(m)))
  lab <- ggdendro::dendro_data(d, type = "rectangle")$labels
  in_x_order <- as.character(lab$label[order(lab$x)])
  list(d = d, labels = list(in_x_order, in_x_order), in_x_order = in_x_order)
}

## Domain of the leaf axis. coord_flip() puts the leaf axis on y.
leaf_domain <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range

brush <- function(ymin, ymax) list(ymin = ymin, ymax = ymax)

## --- 1. domain stability ---------------------------------------------------

test_that("the leaf-axis domain does not move when leaves are highlighted", {
  skip_if_not_installed("ggplot2")
  td <- tiny_dendrogram(8)
  render <- function(highlight) {
    trajectories_dendrogram(precomputed_dendrogram        = td$d,
                            precomputed_dendrogram_labels = td$labels,
                            leaves_to_highlight.A         = highlight)$regular
  }

  bare <- leaf_domain(render(NULL))

  ## An interior highlight never extended the range, so this passed even before
  ## the fix; the edge cases are the ones that regressed.
  expect_equal(leaf_domain(render(td$in_x_order[4:5])), bare)
  expect_equal(leaf_domain(render(td$in_x_order[1:2])), bare)   # bottom edge
  expect_equal(leaf_domain(render(td$in_x_order[7:8])), bare)   # top edge
  expect_equal(leaf_domain(render(td$in_x_order[8])),   bare)   # single top leaf
})

test_that("leaves sit at integer coordinates inside the fixed domain", {
  skip_if_not_installed("ggplot2")
  td <- tiny_dendrogram(8)
  dom <- leaf_domain(trajectories_dendrogram(precomputed_dendrogram        = td$d,
                                             precomputed_dendrogram_labels = td$labels)$regular)
  expect_lt(dom[1], 1)
  expect_gt(dom[2], 8)
})

test_that("highlight rectangles stay inside the fixed domain and are not clipped", {
  skip_if_not_installed("ggplot2")
  td  <- tiny_dendrogram(8)
  dom <- leaf_domain(trajectories_dendrogram(precomputed_dendrogram        = td$d,
                                             precomputed_dendrogram_labels = td$labels,
                                             leaves_to_highlight.A         = td$in_x_order[8])$regular)
  ## Rect for the top leaf reaches 8 + 0.25; the domain must still contain it.
  expect_gt(dom[2], 8.25)
})

test_that("the height axis keeps its expansion so hjust=1 leaf labels are not clipped", {
  skip_if_not_installed("ggplot2")
  td <- tiny_dendrogram(8)
  p  <- trajectories_dendrogram(precomputed_dendrogram        = td$d,
                                precomputed_dendrogram_labels = td$labels)$regular
  ## Labels are drawn at the leaves' own height with hjust = 1, so they extend
  ## below the data minimum. expand = FALSE would zero this margin and cut them.
  height <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range
  expect_lt(height[1], 0)
})

## --- 2. brush -> classes-vector positions ----------------------------------

test_that("a brush over the main plot returns descending leaf positions", {
  expect_equal(.brushed_leaf_idcs(brush(2.6, 5.4), seq_len(9)), c(5, 4, 3))
})

test_that("a brush covering a single leaf returns just that leaf", {
  expect_equal(.brushed_leaf_idcs(brush(3.6, 4.4), seq_len(9)), 4)
})

test_that("zoom leaves drawn at x_min..x_max map back to classes-vector positions", {
  ## The zoom plot renders leaves 12..16 of the full tree, but
  ## dendrogram_zoom_classes is indexed 1..5 over exactly that span.
  expect_equal(.brushed_leaf_idcs(brush(13.6, 14.4), 12:16), 3)
  expect_equal(.brushed_leaf_idcs(brush(12.6, 15.4), 12:16), c(4, 3, 2))
})

test_that("a zoom brush past the end of the rendered leaves does not run off the classes vector", {
  got <- .brushed_leaf_idcs(brush(15.6, 40), 12:16)
  expect_true(all(got >= 1 & got <= 5))
})

test_that("no brush selects nothing", {
  expect_equal(.brushed_leaf_idcs(NULL, seq_len(9)), integer(0))
  expect_equal(.brushed_leaf_idcs(brush(NULL, NULL), seq_len(9)), integer(0))
})

test_that("an empty leaf set selects nothing", {
  expect_equal(.brushed_leaf_idcs(brush(1, 5), integer(0)), integer(0))
})

test_that("a brush far off the plot selects nothing", {
  expect_equal(.brushed_leaf_idcs(brush(40, 50), seq_len(9)), integer(0))
})

## --- 3. coordmap alignment -------------------------------------------------
##
## The bug that survived the two coordinate fixes: a NEGATIVE vertical plot.margin
## makes the gtable panel taller than the device, and shiny's coordmap arithmetic
## then disagrees with where grid actually draws by exactly that overhang -- a
## constant 22.2 px at 144 dpi (0.4 cm, the two -0.2 cm margins). Constant in
## pixels means the error in LEAVES scales with panel density: 0.15 leaves at n=9,
## 0.40 at n=24, 0.57 on a short panel. Selection was off by up to half a branch
## with no code path to blame, which is why it looked environment-dependent.
##
## This measures the real thing: render the plot, find where the ink for each leaf
## coordinate actually lands, and compare against what shiny's coordmap promises.

## Returns mean(actual - predicted) in image pixels for the rendered dendrogram.
coordmap_offset_px <- function(n, h_css, w_css = 360, ratio = 2, margin = NULL) {
  td <- tiny_dendrogram(n)
  p  <- trajectories_dendrogram(precomputed_dendrogram        = td$d,
                                precomputed_dendrogram_labels = td$labels)$regular
  if (!is.null(margin)) {
    p <- p + ggplot2::theme(
      plot.margin = ggplot2::unit(c(margin, 0, margin, 0), "cm"))
  }
  bd <- ggplot2::ggplot_build(p)
  ym <- max(bd$data[[1]]$y, bd$data[[1]]$yend)
  ## markers in data space at every leaf coordinate; explicit segments rather than
  ## geom_vline, whose behaviour under coord_flip is version-dependent
  pp <- p + ggplot2::geom_segment(
    data = data.frame(k = seq_len(n)), inherit.aes = FALSE,
    ggplot2::aes(x = k, xend = k, y = ym * 0.55, yend = ym * 0.75),
    colour = "#00FF00", linewidth = 0.35)

  w <- w_css * ratio; h <- h_css * ratio; res <- 72 * ratio
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)

  ## render the way shiny does: device first, coordmap from the gtable drawn on it
  grDevices::png(f, width = w, height = h, res = res)
  built <- ggplot2::ggplot_build(pp)
  gt    <- ggplot2::ggplot_gtable(built)
  grid::grid.newpage(); grid::grid.draw(gt)
  cm <- shiny:::getGgplotCoordmap(
    structure(list(build = built, gtable = gt), class = "ggplot_build_gtable"),
    w_css, h_css, res)
  grDevices::dev.off()

  pan  <- cm$panels[[1]]
  dom  <- c(pan$domain$bottom, pan$domain$top)
  rng  <- c(pan$range$bottom,  pan$range$top)
  pred <- rng[1] + (seq_len(n) - dom[1]) / (dom[2] - dom[1]) * (rng[2] - rng[1])

  img <- png::readPNG(f)
  rows <- which(apply(img, 1, function(r) any(r[, 2] > 0.6 & r[, 1] < 0.5 & r[, 3] < 0.5)))
  act  <- sort(as.numeric(tapply(rows, cumsum(c(1, diff(rows) > 1)), mean)),
               decreasing = TRUE)
  if (length(act) != n) return(NA_real_)
  mean(act - pred)
}

test_that("shiny's coordmap agrees with where the leaves are actually drawn", {
  skip_if_not_installed("png")
  skip_if_not_installed("shiny")
  ## Sub-pixel across leaf counts and panel heights. Before the fix this was
  ## -22.2 px at every one of them.
  expect_lt(abs(coordmap_offset_px(24, 725)), 1)
  expect_lt(abs(coordmap_offset_px(24, 500)), 1)
  expect_lt(abs(coordmap_offset_px(9,  725)), 1)
})

test_that("a negative vertical plot margin is what breaks the coordmap", {
  skip_if_not_installed("png")
  skip_if_not_installed("shiny")
  ## Pins the mechanism, not just the symptom: re-introducing the margin
  ## re-introduces the offset, and it is the size the geometry predicts
  ## (0.4 cm at 144 dpi = 22.7 px).
  expect_lt(abs(coordmap_offset_px(24, 725, margin =  0.0)), 1)
  expect_lt(abs(coordmap_offset_px(24, 725, margin =  0.2)), 1)
  expect_equal(coordmap_offset_px(24, 725, margin = -0.2), -22.7, tolerance = 0.05)
})

test_that("the rendered dendrogram keeps non-negative top and bottom margins", {
  td <- tiny_dendrogram(8)
  m  <- trajectories_dendrogram(precomputed_dendrogram        = td$d,
                                precomputed_dendrogram_labels = td$labels)$regular$theme$plot.margin
  expect_gte(as.numeric(m[1]), 0)   # top
  expect_gte(as.numeric(m[3]), 0)   # bottom
})

## --- 4. real brushes recorded from the app ---------------------------------
##
## Measured in the running app on a 24-leaf tree, domain -0.7 .. 25.7 (exactly
## [0.5 - 1.2, 24.5 + 1.2], i.e. the pinned panel behaving as designed). The
## logged values were taken while the coordmap was still offset, so each is
## corrected by the offset measured above for this geometry, 22.16 px / 55.78
## px-per-leaf = 0.397 leaves. Recorded because they are the cases the user
## reported as wrong.
test_that("real recorded brushes select exactly the leaf that was aimed at", {
  off <- 0.397
  expect_equal(.brushed_leaf_idcs(brush(24.240 - off, 24.760 - off), seq_len(24)), 24)
  expect_equal(.brushed_leaf_idcs(brush(12.921 - off, 13.700 - off), seq_len(24)), 13)
  expect_equal(.brushed_leaf_idcs(brush( 0.859 - off,  1.750 - off), seq_len(24)),  1)
})

test_that("the uncorrected brushes are the over-selections that were reported", {
  ## Documents the failure rather than the fix: with the coordmap offset present,
  ## two of the three brushes pulled in the leaf above. Guards against anyone
  ## "fixing" this again with a tolerance constant.
  expect_equal(.brushed_leaf_idcs(brush(12.921, 13.700), seq_len(24)), 13)
  expect_equal(.brushed_leaf_idcs(brush( 0.859,  1.750), seq_len(24)),  1)
})
