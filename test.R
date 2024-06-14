library(rvinecopulib)

u <- rbicop(1000, "clayton", 0, 3)

quadtree:::test(u, rbind(c(0.5, 0.8), c(0.5, 1)), n_samples = 10, depth = 8)


microbenchmark::microbenchmark(
  quadtree:::test(u, rbind(c(0.5, 0.8), c(0.5, 1)), n_samples = 100, depth = 8),
  times = 20
)
