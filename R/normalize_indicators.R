#' @title Normalize indicators
#'
#' @description Normalize the values of the indicators according to Shen at al. (2013).
#'
#' @usage normalize_indicators(data_indicators,
#'                      cost = NULL)
#'
#' @param data_indicators A data frame with the values of indicators of each DMU by rows. The first column
#'  contains the names of the DMUs.
#' @param cost An array with the positions of cost indicators (1 for the first indicator, 2 for the
#'  second one, etc.). The rest are supposed to be benefit indicators. By default, there are no cost
#'  indicators.
#'
#' @returns A data frame with the values of normalized indicators of each DMU by rows.
#'
#' @author
#' \strong{Vicente Bolós} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Rafael Benítez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @references
#' Yongjun Shen; Elke Hermans; Tom Brijs; Geert Wets (2013). "Data Envelopment Analysis
#' for Composite Indicators: A Multiple Layer Model", Social Indicators Research 114, 739–756.
#' \doi{10.1007/s11205-012-0171-0}
#'
#' @examples
#' # Data indicators in Shen et al. (2013) are, in fact, not normalized.
#'
#' data("Shen2013")
#' Shen2013_normalized <- normalize_indicators(Shen2013)
#'
#' @seealso \code{\link{model_multilayer}}, \code{\link{model_basicCI}}
#'
#' @export


normalize_indicators <-
  function(data_indicators,
           cost = NULL) {

    Y <- as.matrix(data_indicators[, -1])
    nI <- ncol(Y)     # Número de Indicadores

    if (!is.null(cost)) {
      cost <- sort(unique(cost))
      if (!all(cost %in% 1:nI)) {
        stop("Parameter cost must correspond to positions of cost indicators (1 for the first indicator, 2 for the second one, etc.)")
      }
    }
    benefit <- setdiff(1:nI, cost)

    Yb <- Y[, benefit, drop = FALSE]
    Yc <- Y[, cost, drop = FALSE]
    max_cols <- apply(Yb, 2, max)
    min_cols <- apply(Yc, 2, min)

    Y_norm <- matrix(0, nrow = nrow(Y), ncol = nI)
    Y_norm[, benefit] <- sweep(Yb, 2, max_cols, FUN = "/")
    Y_norm[, cost] <- sweep(Yc, 2, min_cols, FUN = "/")
    Y_norm[, cost] <- 1 / Y_norm[, cost]

    data_indicators_norm <- cbind(data_indicators[, 1], Y_norm)
    colnames(data_indicators_norm) <- colnames(data_indicators)

    return(data_indicators_norm)

  }
