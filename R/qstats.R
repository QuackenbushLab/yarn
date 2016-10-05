#' Compute quantile statistics
#'
#' This function was directly borrowed from github user kokrah.
#'
#' @param exprs for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param window window size for running median as a fraction on the number of rows of exprs
#'
#' @importFrom stats runmed
#' @importFrom stats model.matrix
#'
#' @return list of statistics
#'
#' @source \href{https://raw.githubusercontent.com/kokrah/qsmooth/master/R/qstats.r}{Kwame Okrah's qsmooth R package}
#' Compute quantile statistics
#'
qstats <- function(exprs, groups, window) {
  # Compute sample quantiles
  Q <- apply(exprs, 2, sort)
  # Compute quantile reference
  Qref <- rowMeans(Q)
  # Compute SST
  SST <- rowSums((Q - Qref)^2)
  # Compute SSB
  f <- factor(as.character(groups))
  X <- model.matrix(~0 + f)
  QBETAS <- t(solve(t(X) %*% X) %*% t(X) %*% t(Q))
  Qhat <- QBETAS %*% t(X)
  SSB <- rowSums((Qhat - Qref)^2)
  # Compute weights
  roughWeights <- 1 - SSB/SST
  roughWeights[SST < 1e-06] <- 1
  # Compute smooth weights
  k <- floor(window * nrow(Q))
  if (k%%2 == 0)
    k <- k + 1
  smoothWeights <- runmed(roughWeights, k = k, endrule = "constant")
  list(Q = Q, Qref = Qref, Qhat = Qhat, QBETAS = QBETAS, SST = SST,
       SSB = SSB, SSE = SST - SSB, roughWeights = roughWeights,
       smoothWeights = smoothWeights)
}
