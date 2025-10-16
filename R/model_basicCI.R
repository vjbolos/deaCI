#' @title Basic DEA-CI model
#'
#' @description Basic DEA model for composite indicators, according to Shen at al. (2013).
#'
#' @usage model_basicCI(data_indicators,
#'               ubounds = NULL,
#'               dmu_eval = NULL,
#'               dmu_ref = NULL,
#'               returnlp = FALSE)
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
#' result <- model_basicCI(data_indicators = Shen2013)
#' CI(result)
#'
#' @seealso \code{\link{model_basicCI_super}}, \code{\link{model_multilayer}}, \code{\link{CI}}
#'
#' @import lpSolve
#'
#' @export

model_basicCI <-
  function(data_indicators,
           ubounds = NULL,
           dmu_eval = NULL,
           dmu_ref = NULL,
           returnlp = FALSE) {

    Y <- as.matrix(data_indicators[, -1])
    nd <- nrow(Y)     # Número de DMUs
    nI <- ncol(Y)     # Número de Indicadores
    dmunames <- data_indicators[, 1]
    Inames <- names(data_indicators[, -1])

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

    # --- Modelo Basic DEA-CI ---

    obj <- "max"

    # Constraints of 1st block
    f.con.1 <- Y[dmu_ref, ]
    f.dir.1 <- rep("<=", ndr)
    f.rhs.1 <- rep(1, ndr)
    if (is.null(ubounds)) {
      # All Constraints
      f.con <- f.con.1
      f.dir <- f.dir.1
      f.rhs <- f.rhs.1
    } else {
      if ((length(ubounds) != 2) || (ubounds[1] > (1 / nI)) || (ubounds[2] < (1 / nI)) || (ubounds[2] > 1) || (ubounds[1] < 0)) {
        stop("Parameter ubounds incorrect.")
      }
    }

    DMU <- vector("list", nde)
    names(DMU) <- dmunames[dmu_eval]
    aux <- rep(0, nI)

    for (dmu_idx in 1:nde) {

      ii <- dmu_eval[dmu_idx]

      f.obj <- Y[ii, ]   # Función objetivo

      # Constraints of 3rd Block (Cotas para u*y_0)

      if (!is.null(ubounds)) {
        f.con.umin <- NULL
        f.con.umax <- NULL
        aux21 <- rep(ubounds[1], nI) * Y[ii, ]
        aux22 <- rep(ubounds[2], nI) * Y[ii, ]
        for (f in 1:nI) {
          aux1 <- aux
          aux1[f] <- Y[ii, f]
          f.con.umin <- rbind(f.con.umin, aux1 - aux21)
          f.con.umax <- rbind(f.con.umax, aux1 - aux22)
        }
        f.dir.umin <- rep(">=", nI)
        f.dir.umax <- rep("<=", nI)
        f.rhs.umin <- rep(0, nI)
        f.rhs.umax <- f.rhs.umin

        # All Constraints
        f.con <- rbind(f.con.1, f.con.umin, f.con.umax)
        f.dir <- c(f.dir.1, f.dir.umin, f.dir.umax)
        f.rhs <- c(f.rhs.1, f.rhs.umin, f.rhs.umax)
      }

      if (returnlp) {
        u <- rep(0,nI)
        names(u) <- Inames
        var = list(u = u)
        DMU[[dmu_idx]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con,
                               const.dir = f.dir, const.rhs = f.rhs, var = var)
      } else {

        res <- lp(obj, f.obj, f.con, f.dir, f.rhs)

        if (res$status == 0) {
          CI <- res$objval
          u <- res$solution
          names(u) <- Inames
        } else {
          CI <- NA
          u <- NA
        }
        DMU[[dmu_idx]] <- list(CI = CI, u = u)
      }

    }

    modelOutput <- list(modelname = "basicCI",
                        DMU = DMU,
                        data_indicators = data_indicators,
                        ubounds = ubounds,
                        dmu_eval = dmu_eval,
                        dmu_ref = dmu_ref
    )

    return(modelOutput)

  }
