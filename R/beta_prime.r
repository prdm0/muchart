#' Calculate the sample mean of a Beta Prime distribution
#'
#' @importFrom extraDistr dbetapr
#' @param x A numeric vector of observations
#' @param n The number of observations in each sample
#' @param mu The mean parameter of the Beta Prime distribution
#' @param phi The shape parameter of the Beta Prime distribution
#'
#' @return The sample mean of the Beta Prime distribution
#'
#' @examples
#' x <- c(0.2, 0.3, 0.4)
#' d_sample_mean_bp(x, n = 3, mu = 1, phi = 2)
#'
#' @export
d_sample_mean_bp <- function(x, n, mu, phi) {
  alpha <- mu * (phi + 1)
  beta <- phi + 2

  lambda <- n * alpha * (alpha + beta^2 - 2 * beta + n * alpha * beta - 2 * n * alpha + 1) / ((beta - 1) * (alpha + beta - 1))
  delta <- (2 * alpha + beta^2 - beta + n * alpha * beta - 2 * n * alpha) / (alpha + beta - 1)

  n * extraDistr::dbetapr(x = n * x, shape1 = lambda, shape2 = delta)
}

#' Generate random samples from a Beta Prime distribution
#'
#' @importFrom extraDistr rbetapr
#'
#' @param lots The number of samples to generate
#' @param n The number of observations in each sample
#' @param mu The mean parameter of the Beta Prime distribution
#' @param phi The shape parameter of the Beta Prime distribution
#' @param ... Additional arguments to be passed to extraDistr::rbetapr
#'
#' @return A matrix of random samples from the Beta Prime distribution
#'
#' @examples
#' r_bp(lots = 10, n = 5, mu = 1, phi = 2)
#'
#' @export
r_bp <- function(lots, n, mu, phi, ...) {
  observations <- function() {
    alpha <- mu * (phi + 1)
    beta <- phi + 2
    extraDistr::rbetapr(n = n * lots, shape1 = alpha, shape2 = beta, ...)
  }
  m <- matrix(data = observations(), nrow = lots, ncol = n, byrow = TRUE)
  colnames(m) <- paste("n", 1L:ncol(m), sep = "_")
  rownames(m) <- paste("sample", 1L:nrow(m), sep = "_")
  return(m)
}

#' @importFrom stats sd
#' @importFrom numDeriv grad
#' @importFrom lbfgs lbfgs
mle_bp <-
  function(data,
           start = c(1, 1),
           epsilon = 1e-6,
           linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
           ...) {
    x <- data
    n <- ncol(x)
    x <- apply(X = x, FUN = \(i) mean(i), MARGIN = 1L)

    log_like_bp <- function(x, n, mu, phi) {
      -sum(log(d_sample_mean_bp(x = x, n = n, mu = mu, phi = phi)))
    }

    grad_bp <- function(x, n, mu, phi) {
      numDeriv::grad(
        func = \(par) log_like_bp(x = x, n = n, mu = par[1L], phi = par[2L]),
        x = c(mu, phi)
      )
    }

    lbfgs::lbfgs(
      call_eval = \(par) log_like_bp(x = x, n = n, mu = par[1L], phi = par[2L]),
      call_grad = \(par) grad_bp(x = x, n = n, mu = par[1L], phi = par[2L]),
      vars = start,
      invisible = 1,
      linesearch_algorithm = linesearch_algorithm,
      epsilon = epsilon,
      ...
    )
  }

#' @importFrom extraDistr qbetapr
q_bp <- function(p, n, mu, phi) {
  alpha1 <- mu * (phi + 1)
  beta <- phi + 2
  lambda <- n * alpha1 * (alpha1 + beta^2 - 2 * beta + n * alpha1 * beta - 2 * n * alpha1 + 1) / ((beta - 1) * (alpha1 + beta - 1))
  delta <- (2 * alpha1 + beta^2 - beta + n * alpha1 * beta - 2 * n * alpha1) / (alpha1 + beta - 1)
  extraDistr::qbetapr(p = p, shape1 = lambda, shape2 = delta) / n
}

limits_bp <-
  function(data = NULL,
           alpha = 0.0027,
           mu = NULL,
           phi = NULL,
           ...) {
    if (!is.matrix(data)) {
      stop("The object data must be a matrix!")
    }

    if (is.null(mu) || is.null(phi)) {
      mle <- mle_bp(data = data, ...)
      mu <- mle$par[1L]
      phi <- mle$par[2L]
      convergence <- mle$convergence
      value <- mle$value
    } else {
      convergence <- NULL
      value <- NULL
    }

    n <- ncol(data)

    r <- q_bp(p = c(alpha / 2, 1 - alpha / 2), n = n, mu = mu, phi = phi)
    return(
      list(
        li = r[1L],
        ls = r[2L],
        mu_hat = mu,
        phi_hat = phi,
        convergence = convergence,
        value = value
      )
    )
  }

#' Calculate statistical measures for a Beta Prime control chart
#'
#' @param data A matrix of observed data
#' @param alpha The significance level for the control limits
#' @param mu The mean parameter of the Beta Prime distribution (optional)
#' @param phi The shape parameter of the Beta Prime distribution (optional)
#' @param ... Additional arguments to be passed to limits_bp
#'
#' @return A list containing the estimated parameters, control limits, and statistical measures
#'
#' @examples
#' set.seed(0)
#' dados <- r_bp(lots = 10, n = 5, mu = 1, phi = 2)
#' stats_bp(data = dados, alpha = 0.0027)
#'
#' @export
stats_bp <- function(data, alpha = 0.0027, mu = NULL, phi = NULL, ...) {
  if (!is.matrix(data)) {
    stop("The object data must be a matrix!")
  }
  n <- ncol(data)
  vec_sum <- apply(X = data, FUN = \(i) mean(i), MARGIN = 1L)
  limits <- limits_bp(data = data, alpha = alpha, mu = mu, phi = phi, ...)
  outside <- sum(vec_sum < limits$li) + sum(vec_sum > limits$ls)
  alpha_hat <- outside / length(vec_sum)
  list(
    mu_hat = limits$mu_hat,
    phi_hat = limits$phi_hat,
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
