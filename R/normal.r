#' Calculate the density of the sample mean for a normal distribution
#'
#' @importFrom stats dnorm
#'
#' @param x The value at which to evaluate the density
#' @param n The sample size
#' @param mu The mean of the distribution
#' @param sd The standard deviation of the distribution
#'
#' @return The density of the sample mean
#'
#' @examples
#' d_sample_mean_n(x = 2, n = 10, mu = 0, sd = 1)
#'
#' @export
d_sample_mean_n <- function(x, n, mu, sd) {
  dnorm(x, mean = mu, sd = sqrt(sd^2 / n))
}

#' Generate random samples from a normal distribution
#'
#' @importFrom stats rnorm
#' @param lots The number of lots to generate
#' @param n The sample size
#' @param mu The mean of the distribution
#' @param sd The standard deviation of the distribution
#'
#' @return A matrix of random samples
#'
#' @examples
#' r_n(lots = 100, n = 10, mu = 0, sd = 1)
#'
#' @export
r_n <- function(lots, n, mu, sd) {
  observations <- function() {
    rnorm(n = n * lots, mean = mu, sd = sd)
  }

  m <- matrix(data = observations(), nrow = lots, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("lot", 1L:nrow(m), sep = "_")
  return(m)
}

#' @importFrom stats sd
#' @importFrom numDeriv grad
#' @importFrom lbfgs lbfgs
mle_n <- function(data, start = c(1, 1), epsilon = 1e-6, linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO", ...) {
  x <- data
  n <- ncol(x)
  x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

  log_like_n <- function(x, n, mu, sd) {
    -sum(log(d_sample_mean_n(x = x, n = n, mu = mu, sd = sd)))
  }

  grad_n <- function(x, n, mu, sd) {
    numDeriv::grad(
      func = \(par) log_like_n(x = x, n = n, mu = par[1L], sd = par[2L]),
      x = c(mu, sd)
    )
  }

  lbfgs::lbfgs(
    call_eval = \(par) log_like_n(x = x, n = n, mu = par[1L], sd = par[2L]),
    call_grad = \(par) grad_n(x = x, n = n, mu = par[1L], sd = par[2L]),
    vars = start,
    invisible = 1L,
    linesearch_algorithm = linesearch_algorithm,
    epsilon = epsilon,
    ...
  )
}

#' @importFrom stats qnorm
limits_n <- function(data = NULL, alpha = 0.0027, mu = NULL, sd = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }

  if (is.null(mu) || is.null(sd)) {
    mle <- mle_n(data = data, ...)
    mu <- mle$par[1L]
    sd <- mle$par[2L]
    convergence <- mle$convergence
    value <- mle$value
  } else {
    convergence <- NULL
    value <- NULL
  }

  n <- ncol(data)
  r <- qnorm(p = c(alpha / 2, 1 - alpha / 2), mean = mu, sd = sqrt(sd^2 / n))

  return(
    list(
      li = r[1L],
      ls = r[2L],
      mu_hat = mu,
      sd_hat = sd,
      convergence = convergence,
      value = value
    )
  )
}

#' Calculate statistics for a normal distribution control chart
#'
#' @param data The data matrix
#' @param alpha The significance level
#' @param mu The mean of the distribution (optional)
#' @param sd The standard deviation of the distribution (optional)
#' @param ... Additional arguments to be passed to the optimization algorithm
#'
#' @return A list containing the estimated parameters, control limits, alpha_hat, ARL, MRL, SDRL, convergence status, and value of the objective function
#'
#' @examples
#' set.seed(0)
#' dados <- r_n(lots = 100e3L, n = 25, mu = 1, sd = 1.7)
#' stats_n(data = dados, alpha = 0.0027)
#'
#' @export
stats_n <- function(data, alpha = 0.0027, mu = NULL, sd = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limits <- limits_n(data = data, alpha = alpha, mu = mu, sd = sd, ...)
  outside <- sum(vec_sum < limits$li) + sum(vec_sum > limits$ls)
  alpha_hat <- outside / length(vec_sum)
  list(
    mu_hat = limits$mu_hat,
    sd_hat = limits$sd_hat,
    alpha_hat = alpha_hat,
    ARL = 1 / alpha_hat,
    MRL = log(0.5) / log(1 - alpha_hat),
    SDRL = sqrt((1 - alpha_hat) / (alpha_hat^2)),
    li = limits$li,
    ls = limits$ls,
    convergence = limits$convergence,
    value = limits$value
  )
}
