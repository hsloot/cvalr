suniroot <- function(
    f, interval, is_increasing, ...,
    lower = min(interval), upper = max(interval),
    tolerance = .Machine$double.eps^0.25,
    max_iter = 1e3,
    p = 0.6, grid_length = 1e2) {
  interval <- c(lower, upper)
  stopifnot(isTRUE(all(is.finite(interval))))
  normalize <- function(x) { x / sum(x) }

  atoms <- seq(lower, upper, length.out = grid_length)
  weights <-  normalize(rep(1, grid_length))
  qf <- q.r(DiscreteDistribution(supp = atoms, prob = weights))

  iter <- 0
  while (!isTRUE(abs(qf(1 - 0.05/2) - qf(0.05/2)) < tolerance) && iter < max_iter) {
    x <- qf(0.5)
    z <- f(x, ...)
    y <- ((z < 0) * 2 - 1) * (is_increasing * 2 - 1)

    i <- min(which(atoms >= x))
    atoms <- c(atoms, x)
    weights <- c(weights, mean(weights[c(i, i+1)]))
    sort_ind <- order(atoms)
    atoms <- atoms[sort_ind]
    weights <- weights[sort_ind]
    weights[y*(atoms  - x) > 0] <- weights[y*(atoms  - x) > 0] * p
    weights[y*(atoms - x) <= 0] <- weights[y*(atoms - x) <= 0] * (1-p)
    weights <- normalize(weights)
    qf <- q.r(DiscreteDistribution(supp = atoms, prob = weights))
    iter <- iter + 1
  }

  qf(0.5)
}
