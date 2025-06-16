#' Robust M-quantile periodogram
#' This function computes the univariate robust M-quantile periodogram using
#' M-quantile regression.
#' @param series univariate time series
#' @param tau quantile level (0 < tau < 1)
#' @return a numeric vector containing the robust estimates of the spectral density
#' @author Patrick Ferreira Patrocinio

mqper <- function(series, tau) {
  n <- length(series)
  perior <- FFT <- NULL
  g <- n %/% 2
  for (j in 1:g) {
    X1 <- X2 <- NULL
    w <- 2 * pi * j / n
    for (i in 1:n) {
      X1[i] <- cos(w * i)
      X2[i] <- sin(w * i)
    }
    if (j != (n / 2)) {
      MX <- cbind(X1, X2)
      fitrob <- mqlm(MX, series, maxit=100, q=tau, k=1.345)
      FFT[j] <- sqrt(n / (8 * pi)) * complex(real = fitrob$coef[1], imaginary = -fitrob$coef[2])
    }
    else {
      MX <- cbind(X1)
      fitrob <- mqlm(MX, series, maxit=100, q=tau, k=1.345)
      FFT[j] <- sqrt(n / (2 * pi)) * complex(real = fitrob$coef[1], imaginary = -0)
    }
    perior[j] <- Mod(FFT[j])^2
  }
  if ((n %% 2) != 0) {
    spec <- c(perior, rev(perior))
  } else {
    spec <- c(perior, rev(perior))
  }
  w <- 2 * pi * seq.int(g) / n
  return(list(spec = spec, freq = w))
}

