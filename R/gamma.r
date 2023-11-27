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
#'
#' set.seed(0)
#' x <- r_g(lots = 1000, n = 10, mu = 1, k = 1.7)
#'
#' # Estimate for maximum likelihood
#' stats_g(data = x)
#'
#' # Using the true parameters
#' stats_g(data = x, mu = 1, k = 1.7)
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

#' Control chart to monitor sample mean of gamma observations
#'
#' @importFrom ggplot2 ggplot aes margin geom_point geom_hline scale_color_manual labs guides guide_legend theme element_text
#'
#' @importFrom scattermore geom_scattermore
#'
#' @importFrom cowplot ggdraw add_sub
#'
#' @param data The data points (matrix)
#' @param alpha The significance level for control limits
#' @param mu The mean of the distribution
#' @param k The shape parameter
#'
#' @return An ggplot2 object with the control limits plot
#'
#' @examples
#' # data <- r_g(500, 5, 1, 1.7)
#' # chart_g(data, 0.0027, 1, 1.7)
#'
#' @export
chart_g <- function(data, alpha = 0.0027, mu = NULL, k = NULL) {
  limits <- stats_g(data = data, alpha = alpha, mu = mu, k = k)
  li <- limits$li
  ls <- limits$ls
  x <- 1:nrow(data)
  y <- apply(X = data, MARGIN = 1, FUN = mean)

  if (is.null(mu) || is.null(k)) {
    mu <- limits$mu_hat
    k <- limits$k_hat
  }

  data <- data.frame(
    Observation = x,
    y = y,
    li = li,
    ls = ls,
    mu = mu,
    outside = y < li | y > ls
  )

  p <-
    data |>
    ggplot()

  if (nrow(data) > 100e3L) {
    p <- p + geom_scattermore(aes(x = data$Observation, y = y, color = data$outside))
  } else {
    p <- p + geom_point(aes(x = data$Observation, y = y, color = data$outside), size = 2, alpha = 0.7)
  }

  p <- p + geom_hline(aes(yintercept = li), color = "#00740e", size = 2, alpha = 0.7) +
    geom_hline(aes(yintercept = ls), color = "#00740e", size = 2, alpha = 0.7) +
    scale_color_manual(
      values = c("#5555ff", "#fd3b3b"),
      breaks = c(FALSE, TRUE),
      labels = c("Under control", "Out of control")
    ) +
    geom_hline(aes(yintercept = mu), color = "black", size = 2, alpha = 0.7) +
    labs(
      title = "Control chart for sample mean",
      subtitle = "",
      x = bquote(bold("Observations")),
      y = bquote(bold("Sample mean")),
      color = "Observations"
    ) +
    guides(color = guide_legend(title = NULL)) +
    theme(
      plot.title = element_text(size = 15, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(size = 8),
      axis.title.x = element_text(
        size = 10, face = "bold",
        margin = margin(30, 0, 0, 0)
      ),
      axis.title.y = element_text(
        size = 10, face = "bold",
        margin = margin(0, 30, 0, 0)
      ),
      axis.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.position = "top",
      plot.margin = margin(0, 0, 0, 0, "cm") # Adjust the position of the second annotation
    )

  p <-
    add_sub(
      p,
      bquote(hat(mu) == .(limits$mu)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )

  p <-
    add_sub(
      p,
      bquote(hat(alpha) == .(format(round(limits$alpha_hat, 5), nsmall = 2))),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )
  p <-
    add_sub(
      p,
      paste("UCL = ", format(round(limits$ls, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )

  p <-
    add_sub(
      p,
      paste("LCL = ", format(round(limits$li, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )

  p <-
    add_sub(
      p,
      paste("ARL = ", format(round(limits$ARL, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )

  p <-
    add_sub(
      p,
      paste("MLR = ", format(round(limits$MRL, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )

  p <-
    add_sub(
      p,
      paste("SDLR = ", format(round(limits$SDRL, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 8
    )
  ggdraw(p)
}
