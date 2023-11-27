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
#' x <- r_bp(lots = 1000, n = 10, mu = 1, phi = 1.7)
#'
#' # Estimate for maximum likelihood
#' stats_bp(data = x, alpha = 0.0027)
#'
#' # Useing the true parameters
#' stats_bp(data = x, alpha = 0.0027, mu = 1, phi = 1.7)
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

#' Control chart to monitor sample mean of beta prime observations
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
#' @param phi The shape parameter
#'
#' @return An ggplot2 object with the control limits plot
#'
#' @examples
#' # data <- r_bp(500, 5, 1, 1.7)
#' # chart_bp(data, 0.0027, 1, 1.7)
#'
#' @export
chart_bp <- function(data, alpha = 0.0027, mu = NULL, phi = NULL) {
  limits <- stats_ig(data = data, alpha = alpha, mu = mu, phi = phi)
  li <- limits$li
  ls <- limits$ls
  x <- 1:nrow(data)
  y <- apply(X = data, MARGIN = 1, FUN = mean)

  if (is.null(mu) || is.null(phi)) {
    mu <- limits$mu_hat
    phi <- limits$phi_hat
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
