#' Calculate the inverse Gaussian probability density function
#'
#' This function calculates the probability density function (PDF) of the inverse Gaussian distribution.
#'
#' @importFrom statmod dinvgauss
#'
#' @param x The value at which to evaluate the PDF.
#' @param n The number of observations.
#' @param mu The mean parameter of the distribution.
#' @param lambda The shape parameter of the distribution.
#'
#' @return The PDF value at the given x.
#'
#' @examples
#' d_sample_mean_ig(x = 2, n = 10, mu = 1, lambda = 2)
#'
#' @export
d_sample_mean_ig <- function(x, n, mu, lambda) {
  statmod::dinvgauss(x = x, mean = mu, shape = n * lambda)
}

#' Generate random samples from the inverse Gaussian distribution
#'
#' This function generates random samples from the inverse Gaussian distribution.
#'
#' @importFrom statmod rinvgauss
#'
#' @param lots The number of lots to generate.
#' @param n The number of observations per lot.
#' @param mu The mean parameter of the distribution.
#' @param lambda The shape parameter of the distribution.
#' @param lotes Logical indicating whether to return the lots as separate matrices.
#'
#' @return A matrix of random samples from the inverse Gaussian distribution.
#'
#' @examples
#' r_ig(lots = 5, n = 10, mu = 1, lambda = 2)
#'
#' @export
r_ig <- function(lots, n, mu, lambda, lotes = TRUE) {
  observations <- function() {
    statmod::rinvgauss(n = n * lots, mean = mu, shape = lambda)
  }
  m <- matrix(data = observations(), nrow = lots, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("sample", 1L:nrow(m), sep = "_")
  return(m)
}

#' @importFrom numDeriv grad
#' @importFrom lbfgs lbfgs
mle_ig <-
  function(data,
           start = c(1, 1),
           epsilon = 1e-6,
           linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
           ...) {
    x <- data
    n <- ncol(x)
    x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

    log_like_ig <- function(x, n, mu, lambda) {
      -sum(log(d_sample_mean_ig(x = x, n = n, mu = mu, lambda = lambda)))
    }

    grad_ig <- function(x, n, mu, lambda) {
      numDeriv::grad(
        func = \(par) log_like_ig(x = x, n = n, mu = par[1L], lambda = par[2L]),
        x = c(mu, lambda)
      )
    }

    lbfgs::lbfgs(
      call_eval = \(par) log_like_ig(x = x, n = n, mu = par[1L], lambda = par[2L]),
      call_grad = \(par) grad_ig(x = x, n = n, mu = par[1L], lambda = par[2L]),
      vars = start,
      invisible = 1L,
      linesearch_algorithm = linesearch_algorithm,
      epsilon = epsilon,
      ...
    )
  }

#' @importFrom statmod qinvgauss
limits_ig <-
  function(data = NULL,
           alpha = 0.0027,
           mu = NULL,
           lambda = NULL,
           ...) {
    if (!is.matrix(data)) {
      stop("The object data must be a matrix!")
    }

    if (is.null(mu) || is.null(lambda)) {
      mle <- mle_ig(data = data, ...)
      mu <- mle$par[1L]
      lambda <- mle$par[2L]
      convergence <- mle$convergence
      value <- mle$value
    } else {
      convergence <- NULL
      value <- NULL
    }

    n <- ncol(data)

    r <- statmod::qinvgauss(p = c(alpha / 2, 1 - alpha / 2), mean = mu, shape = n * lambda)

    return(list(
      li = r[1L],
      ls = r[2L],
      mu_hat = mu,
      lambda_hat = lambda,
      convergence = convergence,
      value = value
    ))
  }

#' Calculate statistical measures for the inverse Gaussian control chart
#'
#' This function calculates various statistical measures for the inverse Gaussian control chart.
#'
#' @param data The data matrix.
#' @param alpha The significance level.
#' @param mu The mean parameter of the distribution (optional).
#' @param lambda The shape parameter of the distribution (optional).
#' @param ... Additional arguments to be passed to the limits function.
#'
#' @return A list containing the estimated parameters, alpha-hat, ARL, MRL, SDRL, control limits, convergence status, and value of the objective function.
#'
#' @examples
#' data <- r_ig(lots = 10, n = 20, mu = 1, lambda = 2)
#' stats_ig(data = data)
#'
#' @export
stats_ig <- function(data, alpha = 0.0027, mu = NULL, lambda = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limits <- limits_ig(data = data, alpha = alpha, mu = mu, lambda = lambda, ...)
  outside <- sum(vec_sum < limits$li) + sum(vec_sum > limits$ls)
  alpha_hat <- outside / length(vec_sum)
  list(
    mu_hat = limits$mu_hat,
    lambda_hat = limits$lambda_hat,
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
