#' Calculate the density of the sample mean for a normal distribution
#'
#' @importFrom stats dnorm
#'
#' @param x The value at which to evaluate the density
#' @param n The sample size
#' @param mu The mean of the distribution
#' @param sd The standard deviation of the distribution
#' Construct control limits using echarts4r
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
#' x <- r_n(lots = 1000, n = 10, mu = 1, sd = 1.7)
#'
#' # Estimate for maximum likelihood
#' stats_n(data = x, alpha = 0.0027)
#'
#' # Useing the true parameters
#' stats_n(data = x, alpha = 0.0027, mu = 1, sd = 1.7)
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

#' Control chart to monitor sample mean of normal observations
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
#' @param sd The standard deviation of the distribution
#'
#' @return An ggplot2 object with the control limits plot
#'
#' @examples
#' # data <- r_n(500, 5, 1, 1.7)
#' # chart_n(data, 0.0027, 1, 1.7)
#'
#' @export
chart_n <- function(data, alpha = 0.0027, mu = NULL, sd = NULL) {
  limits <- stats_n(data = data, alpha = alpha, mu = mu, sd = sd)
  li <- limits$li
  ls <- limits$ls
  x <- 1:nrow(data)
  y <- apply(X = data, MARGIN = 1, FUN = mean)

  if (is.null(mu) || is.null(sd)) {
    mu <- limits$mu_hat
    sd <- limits$sd_hat
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
    p <- p + geom_point(aes(x = data$Observation, y = y, color = data$outside), size = 5, alpha = 0.7)
  }

  p <- p + geom_hline(aes(yintercept = li), color = "#00740e", size = 4, alpha = 0.7) +
    geom_hline(aes(yintercept = ls), color = "#00740e", size = 4, alpha = 0.7) +
    scale_color_manual(
      values = c("#5555ff", "#fd3b3b"),
      breaks = c(FALSE, TRUE),
      labels = c("Under control", "Out of control")
    ) +
    geom_hline(aes(yintercept = mu), color = "black", size = 4, alpha = 0.7) +
    labs(
      title = "Normal control chart for sample mean",
      subtitle = "",
      x = bquote(bold("Observations")),
      y = bquote(bold("Sample mean")),
      color = "Observations"
    ) +
    guides(color = guide_legend(title = NULL)) +
    theme(
      plot.title = element_text(size = 30, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(size = 20),
      axis.title.x = element_text(
        size = 20, face = "bold",
        margin = margin(30, 0, 0, 0)
      ),
      axis.title.y = element_text(
        size = 20, face = "bold",
        margin = margin(0, 30, 0, 0)
      ),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
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
      size = 15
    )

  p <-
    add_sub(
      p,
      bquote(hat(alpha) == .(format(round(limits$alpha_hat, 5), nsmall = 2))),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 15
    )
  p <-
    add_sub(
      p,
      paste("UCL = ", format(round(limits$ls, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 15
    )

  p <-
    add_sub(
      p,
      paste("LCL = ", format(round(limits$li, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 15
    )

  p <-
    add_sub(
      p,
      paste("ARL = ", format(round(limits$ARL, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 15
    )

  p <-
    add_sub(
      p,
      paste("MLR = ", format(round(limits$MRL, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 15
    )

  p <-
    add_sub(
      p,
      paste("SDLR = ", format(round(limits$SDRL, 5), nsmall = 2)),
      x = 0,
      hjust = 0,
      vjust = 0,
      size = 15
    )
  ggdraw(p)
}
