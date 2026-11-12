# R1 comment 6: how often a histogram of a single NB(mu = 5, r = 5) sample shows two or more peaks.
# Peak: a value, or a run of equal values (plateau), whose count exceeds the counts on both sides
# (outside the support the count is 0). "Strict" counts a peak only at a single value higher than both neighbours.
# "Large": only peaks whose count is at least 1% of the sample.
set.seed(20261109); R <- 10000
peaks <- function(cnt, strict = FALSE) {
  if (strict) { x <- c(0, cnt, 0); i <- 2:(length(x) - 1); return(which(x[i] > x[i - 1] & x[i] > x[i + 1])) }
  r <- rle(cnt); v <- c(0, r$values, 0); i <- 2:(length(v) - 1)
  st <- cumsum(c(1, r$lengths))[seq_along(r$lengths)]
  st[which(v[i] > v[i - 1] & v[i] > v[i + 1])]
}
for (n in c(270, 411, 4691)) {
  res <- t(replicate(R, {
    y <- rnbinom(n, size = 5, mu = 5); cnt <- tabulate(y + 1)
    a <- peaks(cnt); b <- peaks(cnt, strict = TRUE)
    c(any2 = length(a) >= 2, large2 = sum(cnt[a] >= 0.01 * n) >= 2, any2s = length(b) >= 2, large2s = sum(cnt[b] >= 0.01 * n) >= 2)
  }))
  m <- round(100 * colMeans(res), 1)
  cat(sprintf("n = %4d | two or more peaks: %5.1f%% (strict %5.1f%%) | two or more peaks holding >= 1%%: %5.1f%% (strict %5.1f%%)\n",
              n, m["any2"], m["any2s"], m["large2"], m["large2s"]))
}
