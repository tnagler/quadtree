library(rvinecopulib)
library(microbenchmark)

u <- rbicop(10000, "clayton", 0, 3)

a <- quadtree:::test(u, rbind(c(0.1, 0.1), c(0.5, 0.7)), depth = 6)
a2 <- quadtree:::test2(u, rbind(c(0.1, 0.1), c(0.5, 0.7)), depth = 6)
# b <- quadtree:::test_range(u, rbind(c(0.1, 0.1), c(0.5, 0.7)), depth = 6)
# c(0.5, 0.7)
# b2 <- quadtree:::test_range2(u, rbind(c(0.1, 0.1), c(0.5, 0.7)), depth = 6)
# c(0.5, 0.7)

sum(a - a2)

range <- rbind(c(0.01, 0.01), c(0.9, 0.8))

microbenchmark::microbenchmark(
  quadtree:::test(u, range, depth = 8),
  quadtree:::test2(u, range, depth = 8),
  # quadtree:::test_range(u, range, depth = 8),
  # quadtree:::test_range2(u, range, depth = 8),
  times = 5
)
