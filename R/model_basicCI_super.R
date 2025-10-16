#' @title Radial super-efficiency basic DEA-CI model
#'
#' @description Radial super-efficiency basic DEA model for composite indicators,
#' according to Shen at al. (2013). Radial super-efficiency in the sense of Andersen
#' and Petersen (1993).
#'
#' @usage model_basicCI_super(data_indicators,
#'                     ubounds = NULL,
#'                     dmu_eval = NULL,
#'                     dmu_ref = NULL,
#'                     returnlp = FALSE)
#'
#' @param data_indicators A data frame with the values of indicators of each DMU by rows. The first column
#'  contains the names of the DMUs.
#' @param ubounds A vector with lower and upper relative bounds (in proportion to 1) for \code{u * y_0}.
#'  If \code{NULL} (default), there are no bounds.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective
#'  function and constraints).
#'
#' @returns A list of the results for the evaluated DMUs (\code{DMU} component),
#'  along with any other necessary information to replicate the results, such as the name of the model and
#'  parameters \code{data_indicators}, \code{ubounds}, \code{dmu_eval} and \code{dmu_ref}.
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
#' # Replication of results in Shen et al. (2013) using layer_list.
#'
#' data("Shen2013")
#' result <- model_basicCI_super(data_indicators = Shen2013)
#' CI(result)
#'
#' @seealso \code{\link{model_basicCI}}, \code{\link{model_multilayer_super}}, \code{\link{CI}}
#'
#' @export

model_basicCI_super <-
  function(data_indicators,
           ubounds = NULL,
           dmu_eval = NULL,
           dmu_ref = NULL,
           returnlp = FALSE) {

    nd <- nrow(data_indicators)
    dmunames <- data_indicators[, 1]

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
                             ubounds = ubounds,
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
                        ubounds = ubounds,
                        dmu_eval = dmu_eval,
                        dmu_ref = dmu_ref
    )

    return(modelOutput)

  }
