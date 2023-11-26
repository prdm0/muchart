#' Calculate the gamma distribution for a given sample mean
#'
#' This function calculates the probability density function of the gamma distribution
#' for a given sample mean.
#'
#' @importFrom stats dgamma
#'
#' @param x The sample mean
#' @param n The sample size
#' @param mu The mean of the gamma distribution
#' @param k The shape parameter of the gamma distribution
#' @return The probability density function value
#' @examples
#' d_sample_mean_g(x = 2, n = 10, mu = 5, k = 2)
#'
#' @export
d_sample_mean_g <- function(x, n, mu, k) {
  dgamma(x = x, shape = n * k, scale = mu / (n * k))
}

#' Generate random samples from the gamma distribution
#'
#' This function generates random samples from the gamma distribution.
#'
#' @importFrom stats rgamma
#'
#' @param lots The number of lots
#' @param n The sample size
#' @param mu The mean of the gamma distribution
#' @param k The shape parameter of the gamma distribution
#' @param ... Additional arguments to be passed to `rgamma()` from the `stats` package
#' @return A matrix of random samples
#' @examples
#' r_g(lots = 5, n = 10, mu = 5, k = 2)
#' r_g(lots = 10, n = 20, mu = 10, k = 3)
#'
#' @export
r_g <- function(lots, n, mu, k, ...) {
  observations <- function() {
    rgamma(n = n * lots, shape = k, scale = mu / k, ...)
  }
  m <- matrix(data = observations(), nrow = lots, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("sample", 1L:nrow(m), sep = "_")
  return(m)
}

#' @importFrom stats sd
#' @importFrom numDeriv grad
#' @importFrom lbfgs lbfgs
mle_g <-
  function(data,
           start = c(1, 1),
           epsilon = 1e-6,
           linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
           ...) {
    x <- data
    n <- ncol(x)
    x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

    log_like_g <- function(x, n, mu, k) {
      -sum(log(d_sample_mean_g(x = x, n = n, mu = mu, k = k)))
    }

    grad_g <- function(x, n, mu, k) {
      numDeriv::grad(
        func = \(par) log_like_g(x = x, n = n, mu = par[1L], k = par[2L]),
        x = c(mu, k)
      )
    }

    lbfgs::lbfgs(
      call_eval = \(par) log_like_g(x = x, n = n, mu = par[1L], k = par[2L]),
      call_grad = \(par) grad_g(x = x, n = n, mu = par[1L], k = par[2L]),
      vars = start,
      invisible = 1L,
      linesearch_algorithm = linesearch_algorithm,
      epsilon = epsilon,
      ...
    )
  }

#' @importFrom stats qgamma
limits_g <-
  function(data = NULL,
           alpha = 0.0027,
           mu = NULL,
           k = NULL,
           ...) {
    if (!is.matrix(data)) {
      stop("The object data must be a matrix!")
    }

    if (is.null(mu) || is.null(k)) {
      mle <- mle_g(data = data, ...)
      mu <- mle$par[1L]
      k <- mle$par[2L]
      convergence <- mle$convergence
      value <- mle$value
    } else {
      convergence <- NULL
      value <- NULL
    }

    n <- ncol(data)

    r <- qgamma(p = c(alpha / 2, 1 - alpha / 2), shape = n * k, scale = mu / (n * k))
    return(list(
      li = r[1L],
      ls = r[2L],
      mu_hat = mu,
      k_hat = k,
      convergence = convergence,
      value = value
    ))
  }

#' Calculate statistical measures for the gamma distribution
#'
#' This function calculates statistical measures for the gamma distribution
#' based on the given data.
#'
#' @param data The data matrix
#' @param alpha The significance level
#' @param mu The mean of the gamma distribution (optional)
#' @param k The shape parameter of the gamma distribution (optional)
#' @param ... Additional arguments to be passed `lbfgs::lbfgs()`
#' @return A list containing the statistical measures and other information
#' @examples
#' data <- r_g(lots = 5, n = 10, mu = 5, k = 2)
#' stats_g(data = data)
#'
#' @export
stats_g <- function(data, alpha = 0.0027, mu = NULL, k = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limits <- limits_g(data = data, alpha = alpha, mu = mu, k = k, ...)
  outside <- sum(vec_sum < limits$li) + sum(vec_sum > limits$ls)
  alpha_hat <- outside / length(vec_sum)
  list(
    mu_hat = limits$mu_hat,
    k_hat = limits$k_hat,
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
