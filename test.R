library(rvinecopulib)

u <- rbicop(10000, "clayton", 0, 3)

# quadtree:::test(u, rbind(c(0.1, 0.1), c(0.3, 0.3)), n_samples = 1000, depth = 8)

print(u[1, ])
quadtree:::test(u, rbind(c(0.1, 0.1), c(0.2, 0.3)), n_samples = 100, depth = 7)

microbenchmark::microbenchmark(
  quadtree:::test(u, rbind(c(0.1, 0.1), c(0.7, 0.7)), n_samples = 100, depth = 7),
  times = 5
)

