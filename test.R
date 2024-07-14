library(rvinecopulib)

u <- rbicop(1000, "clayton", 0, 3)

# quadtree:::test(u, rbind(c(0.1, 0.1), c(0.2, 0.3)), n_samples = 100, depth = 6)
a <- quadtree:::test_range(u, rbind(c(0.1, 0.1), c(0.5, 0.7)), depth = 6)


range <- rbind(c(0.01, 0.01), c(0.9, 0.8))

microbenchmark::microbenchmark(
  quadtree:::test(u, range, depth = 8),
  quadtree:::test_range(u, range, depth = 8),
  times = 5
)


