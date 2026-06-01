#!/usr/bin/env Rscript
# =============================================================================
# test_miami_downsample.R  --  the Miami plot thins its dense sub-threshold
# point cloud so a genome-wide plot renders in seconds, NOT minutes. This pins
# downsample_miami_points(): it must (1) keep EVERY peak (abs(y) >= keep_y),
# (2) cap the low band at max_low, (3) be deterministic, and (4) leave the
# global RNG state untouched.
#
# Run from the repo root:   Rscript tests/test_miami_downsample.R
# =============================================================================

source("R/plots.R")

fail <- 0L
check <- function(cond, msg) {
  cat(ifelse(isTRUE(cond), "PASS  ", "FAIL  "), msg, "\n", sep = "")
  if (!isTRUE(cond)) fail <<- fail + 1L
}

set.seed(42)
n_low <- 5e5L
n_high <- 1234L
pts <- data.frame(
  chr = "1",
  pos = seq_len(n_low + n_high),
  # low band: |y| < 4 ; peaks: |y| >= 4 (mix of up/down to mirror both traits)
  y = c(runif(n_low, -3.9, 3.9), sample(c(1, -1), n_high, TRUE) * runif(n_high, 4, 30)),
  trait = "trait 1",
  stringsAsFactors = FALSE
)

# A high-water mark we can search for after thinning.
keep_y <- 4
n_peaks_in <- sum(abs(pts$y) >= keep_y)

out <- downsample_miami_points(pts, keep_y = keep_y, max_low = 2e5L, seed = 1L)

check(nrow(out) < nrow(pts), "downsampled frame is smaller than input")
check(sum(abs(out$y) < keep_y) == 2e5L, "low band capped at max_low")
check(sum(abs(out$y) >= keep_y) == n_peaks_in, "every peak (|y|>=keep_y) kept")
check(nrow(out) == 2e5L + n_peaks_in, "total = capped low band + all peaks")

# Determinism: same seed -> identical rows.
out2 <- downsample_miami_points(pts, keep_y = keep_y, max_low = 2e5L, seed = 1L)
check(identical(out, out2), "same seed gives identical result")

# Global RNG untouched: a draw before and after must match a clean draw.
set.seed(99)
ref <- runif(3)
set.seed(99)
invisible(downsample_miami_points(pts, max_low = 2e5L, seed = 7L))
after <- runif(3)
check(identical(ref, after), "global RNG state is restored")

# No-op when already small enough.
small <- pts[seq_len(1000), , drop = FALSE]
check(identical(downsample_miami_points(small), small), "no-op below max_low")

# Empty input is returned untouched.
empty <- pts[0, , drop = FALSE]
check(identical(downsample_miami_points(empty), empty), "empty input untouched")

if (fail > 0L) {
  cat(sprintf("\n%d check(s) FAILED.\n", fail))
  quit(status = 1L)
}
cat("\nAll miami-downsample checks passed.\n")
