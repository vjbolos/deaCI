#' @title Radial super-efficiency basic DEA-CI model
#'
#' @description Radial super-efficiency basic DEA model for composite indicators,
#' according to Shen at al. (2013). Radial super-efficiency in the sense of Andersen
#' and Petersen (1993).
#'
#' @usage model_basicCI_super(data_indicators,
#'                     sharebounds = NULL,
#'                     ubounds = NULL,
#'                     ubounds_rel = NULL,
#'                     dmu_eval = NULL,
#'                     dmu_ref = NULL,
#'                     returnlp = FALSE)
#'
#' @param data_indicators A data frame with the values of indicators of each DMU by rows. The names of
#'  the DMUs are the names of the rows. Optionally, the first column could contain the names of the DMUs.
#' @param sharebounds It can be a vector of length \code{2} or a matrix of \code{2} rows and a number of
#'  columns equal to the number of indicators. In the first case, it contains
#'  the lower and upper bounds (in proportion to 1) of the shares for all indicators.
#'  The shares are the proportions that each indicator contributes to the efficiency score,
#'  according to Shen et al. (2013). In the second case, it contains the lower and upper bounds
#'  of the shares for each indicator(minimums in the first row and maximums
#'  in the second row).
#'  If \code{sharebounds} is \code{NULL} (default), it is constructed automatically.
#'  For example, if there are 3 categories in the top layer, it takes the value \code{c(0.1, 0.5)},
#'  as in Shen et al. (2013).
#' @param ubounds It can be a vector of length \code{2} or a matrix of \code{2} rows and a number of
#'  columns equal to the number of indicators. In the first case, it contains
#'  the lower and upper bounds of the weights u for all indicators.
#'  In the second case, it contains the lower and upper bounds
#'  of the weights u for each indicator (minimums in the first row and maximums
#'  in the second row). Note that the weights u do not add up to 1.
#' @param ubounds_rel It can be a vector of length \code{2} or a matrix of \code{2} rows and a number of
#'  columns equal to the number of indicators. In the first case, it contains
#'  the lower and upper bounds of the weights u (relativized in proportion to 1) for all indicators.
#'  In the second case, it contains the lower and upper bounds of the weights u (relativized
#'  in proportion to 1) for each indicator (minimums in the first row and maximums
#'  in the second row). Note that the relativized weights sum to 1.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective
#'  function and constraints).
#'
#' @returns A list of the results for the evaluated DMUs (\code{DMU} component),
#'  along with any other necessary information to replicate the results, such as the name of the model and
#'  parameters \code{data_indicators}, \code{sharebounds}, \code{ubounds}, \code{ubounds_rel},
#'  \code{dmu_eval} and \code{dmu_ref}.
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
#' # Replication of results in Shen et al. (2013).
#'
#' data("Shen2013")
#' sharebounds <- c(0, 1)
#' result <- model_basicCI_super(data_indicators = Shen2013,
#'                               sharebounds = sharebounds)
#' CI(result)
#'
#' @seealso \code{\link{model_basicCI}}, \code{\link{model_multilayer_super}}, \code{\link{CI}}
#'
#' @export

model_basicCI_super <-
  function(data_indicators,
           sharebounds = NULL,
           ubounds = NULL,
           ubounds_rel = NULL,
           dmu_eval = NULL,
           dmu_ref = NULL,
           returnlp = FALSE) {

    if (class(data_indicators[[1]]) == "character") {
      rownames(data_indicators) <- data_indicators[[1]]
      data_indicators <- data_indicators[, -1]
    }

    nd <- nrow(data_indicators)
    dmunames <- rownames(data_indicators)

    if (is.null(dmu_eval)) {
      dmu_eval <- 1:nd
    } else if (!all(dmu_eval %in% (1:nd))) {
      stop("Invalid set of DMUs to be evaluated (dmu_eval).")
    }
    names(dmu_eval) <- dmunames[dmu_eval]
    nde <- length(dmu_eval)

    if (is.null(dmu_ref)) {
      dmu_ref <- 1:nd
    } else if (!all(dmu_ref %in% (1:nd))) {
      stop("Invalid set of reference DMUs (dmu_ref).")
    }
    names(dmu_ref) <- dmunames[dmu_ref]
    ndr <- length(dmu_ref)

    model_modelname <- "model_basicCI"

    DMU <- vector(mode = "list", length = nde)
    names(DMU) <- dmunames[dmu_eval]

    for (i in 1:nde) {

      ii <- dmu_eval[i]

      deasol <- do.call(model_modelname,
                        list(data_indicators = data_indicators,
                             sharebounds = sharebounds,
                             ubounds = ubounds,
                             ubounds_rel = ubounds_rel,
                             dmu_eval = ii,
                             dmu_ref = dmu_ref[dmu_ref != ii],
                             returnlp = returnlp
                        )
      )

      DMU[[i]] <- deasol$DMU[[1]]

    }

    modelOutput <- list(modelname = "basicCI_super",
                        DMU = DMU,
                        data_indicators = data_indicators,
                        sharebounds = deasol$sharebounds,
                        ubounds = deasol$ubounds,
                        ubounds_rel = deasol$ubounds_rel,
                        dmu_eval = dmu_eval,
                        dmu_ref = dmu_ref
    )

    return(modelOutput)

  }
