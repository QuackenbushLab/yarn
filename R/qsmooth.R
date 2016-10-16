#' Quantile shrinkage normalization
#'
#' This function was modified from github user kokrah.
#'
#' @param obj for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param norm.factors scaling normalization factors
#' @param plot plot weights? (default=FALSE)
#' @param window window size for running median (a fraction of the number of rows of exprs)
#' @param log Whether or not the data should be log transformed before normalization, TRUE = YES.
#'
#' @importFrom stats ave
#' @importFrom graphics par
#' @importFrom graphics abline
#'
#' @return Normalized expression
#'
#' @source \href{https://raw.githubusercontent.com/kokrah/qsmooth/master/R/qsmooth.r}{Kwame Okrah's qsmooth R package}
#' @examples
#' data(skin)
#' head(yarn:::qsmooth(skin,groups=pData(skin)$SMTSD))
#'
qsmooth <- function(obj, groups, norm.factors = NULL, plot = FALSE,
                    window = 0.05,log=TRUE) {
  stopifnot(class(obj)=="ExpressionSet")
  if(log==TRUE){
    exprs <- log2(exprs(obj)+1)
  } else {
    exprs <- exprs(obj)
  }
  # Stop if exprs contains any NA
  if (any(is.na(exprs)))
    stop("exprs contains NAs (K.Okrah)")
  # Scale normalization step
  if (is.null(norm.factors)) {
    dat <- exprs
  } else {
    dat <- t(t(exprs) - norm.factors)
  }
  # Compute quantile stats
  qs <- qstats(dat, groups, window = window)
  Qref <- qs$Qref
  Qhat <- qs$Qhat
  w <- qs$smoothWeights
  # Weighted quantiles
  normExprs <- w * Qref + (1 - w) * Qhat
  # Re-order normExprs by rank of exprs (columnwise)
  for (i in 1:ncol(normExprs)) {
    # Grab ref. i
    ref <- normExprs[, i]
    # Grab exprs column i
    x <- exprs[, i]
    # Grab ranks of x (using min rank for ties)
    rmin <- rank(x, ties.method = "min")
    # If x has rank ties then average the values of ref at those
    # ranks
    dups <- duplicated(rmin)
    if (any(dups)) {
      # Grab ranks of x (using random ranks for ties) (needed to
      # uniquely identify the indices of tied ranks)
      rrand <- rank(x, ties.method = "random")
      # Grab tied ranks
      tied.ranks <- unique(rmin[dups])
      for (k in tied.ranks) {
        sel <- rrand[rmin == k]  # Select the indices of tied ranks
        ref[sel] <- ave(ref[sel])
      }
    }
    # Re-order ref and replace in normExprs
    normExprs[, i] <- ref[rmin]
  }
  # Plot weights
  if (plot) {
    oldpar <- par(mar = c(4, 4, 1.5, 0.5))
    lq <- length(Qref)
    u <- (1:lq - 0.5)/lq
    if (length(u) > 10000) {
      # do not plot more than 10000 points
      sel <- sample(1:lq, 10000)
      plot(u[sel], w[sel], pch = ".", main = "qsmooth weights",
           xlab = " quantiles", ylab = "Weight", ylim = c(0,
                                                          1))
    } else {
      plot(u, w, pch = ".", main = "qsmooth weights", xlab = "quantiles",
           ylab = "Weight", ylim = c(0, 1))
    }
    abline(h = 0.5, v = 0.5, col = "red", lty = 2)
    par(oldpar)
  }
  rownames(normExprs) <- rownames(exprs)
  colnames(normExprs) <- colnames(exprs)
  normExprs
}
