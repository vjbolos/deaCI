#' @title Basic DEA-CI model
#'
#' @description Basic DEA model for composite indicators, according to Shen at al. (2013).
#'
#' @usage model_basicCI(data_indicators,
#'               sharebounds = NULL,
#'               ubounds = NULL,
#'               ubounds_rel = NULL,
#'               dmu_eval = NULL,
#'               dmu_ref = NULL,
#'               returnlp = FALSE)
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
#' result <- model_basicCI(data_indicators = Shen2013,
#'                         sharebounds = sharebounds)
#' CI(result)
#'
#' @seealso \code{\link{model_basicCI_super}}, \code{\link{model_multilayer}}, \code{\link{CI}}
#'
#' @import lpSolve
#'
#' @export

model_basicCI <-
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

    Y <- as.matrix(data_indicators)
    nd <- nrow(Y)     # Número de DMUs
    nI <- ncol(Y)     # Número de Indicadores
    dmunames <- rownames(data_indicators)
    Inames <- colnames(data_indicators)

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

    # --- RESTRICCIONES DE PESOS ---

    if (!is.null(ubounds)) {

      sharebounds <- NULL
      ubounds_rel <- NULL
      if (length(ubounds) == 2) {
        ubounds <- matrix(c(ubounds[1], ubounds[2]), nrow = 2, ncol = nI)
      } else if ((nrow(ubounds) != 2) || (ncol(ubounds) != nI)) {
        stop("Parameter ubounds incorrect.")
      }
      rownames(ubounds) <- c("min", "max")
      colnames(ubounds) <- Inames

    } else if (!is.null(ubounds_rel)) {

      sharebounds <- NULL
      if (length(ubounds_rel) == 2) {
        ubounds_rel <- matrix(c(ubounds_rel[1], ubounds_rel[2]), nrow = 2, ncol = nI)
      } else if ((nrow(ubounds_rel) != 2) || (ncol(ubounds_rel) != nI)) {
        stop("Parameter ubounds_rel incorrect.")
      }
      rownames(ubounds_rel) <- c("min", "max")
      colnames(ubounds_rel) <- Inames

    } else {

      if (is.null(sharebounds)) {
        if (nI > 2) {
          shareboundsmax <- 1 / (nI - 1)
          shareboundsmin <- shareboundsmax / 5
        } else {
          shareboundsmin <- 0.2
          shareboundsmax <- 0.8
        }
        sharebounds <- matrix(c(shareboundsmin, shareboundsmax), nrow = 2, ncol = nI)
      } else if (length(sharebounds) == 2) {
        sharebounds <- matrix(c(sharebounds[1], sharebounds[2]), nrow = 2, ncol = nI)
      } else if ((nrow(sharebounds) != 2) || (ncol(sharebounds) != nI)) {
        stop("Parameter sharebounds incorrect.")
      }
      rownames(sharebounds) <- c("min", "max")
      colnames(sharebounds) <- Inames

    }

    # --- Modelo Basic DEA-CI ---

    obj <- "max"

    # Constraints of 1st block
    f.con.1 <- Y[dmu_ref, ]
    f.dir.1 <- rep("<=", ndr)
    f.rhs.1 <- rep(1, ndr)

    DMU <- vector("list", nde)
    names(DMU) <- dmunames[dmu_eval]
    aux <- rep(0, nI)

    for (dmu_idx in 1:nde) {

      ii <- dmu_eval[dmu_idx]

      f.obj <- Y[ii, ]   # Función objetivo

      # Constraints of 3rd Block (Shares: cotas para u_hat*y_0)

      f.con.umin <- NULL
      f.con.umax <- NULL
      f.dir.umin <- rep(">=", nI)
      f.dir.umax <- rep("<=", nI)

      if (!is.null(ubounds)) {

        f.rhs.umin <- NULL
        f.rhs.umax <- NULL
        for (f in 1:nI) {
          aux1 <- aux
          aux1[f] <- 1
          f.con.umin <- rbind(f.con.umin, aux1)
          f.rhs.umin <- c(f.rhs.umin, ubounds[1, f])
          f.rhs.umax <- c(f.rhs.umax, ubounds[2, f])
        }
        f.con.umax <- f.con.umin

      } else if (!is.null(ubounds_rel)) {

        for (f in 1:nI) {
          aux21 <- rep(ubounds_rel[1, f], nI)
          aux22 <- rep(ubounds_rel[2, f], nI)
          aux1 <- aux
          aux1[f] <- 1
          f.con.umin <- rbind(f.con.umin, aux1 - aux21)
          f.con.umax <- rbind(f.con.umax, aux1 - aux22)
        }
        f.rhs.umin <- rep(0, nI)
        f.rhs.umax <- f.rhs.umin

      } else {

        for (f in 1:nI) {
          aux21 <- rep(sharebounds[1, f], nI) * Y[ii, ]
          aux22 <- rep(sharebounds[2, f], nI) * Y[ii, ]
          aux1 <- aux
          aux1[f] <- Y[ii, f]
          f.con.umin <- rbind(f.con.umin, aux1 - aux21)
          f.con.umax <- rbind(f.con.umax, aux1 - aux22)
        }
        f.rhs.umin <- rep(0, nI)
        f.rhs.umax <- f.rhs.umin

      }

      # All Constraints
      f.con <- rbind(f.con.1, f.con.umin, f.con.umax)
      f.dir <- c(f.dir.1, f.dir.umin, f.dir.umax)
      f.rhs <- c(f.rhs.1, f.rhs.umin, f.rhs.umax)


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
                        sharebounds = sharebounds,
                        ubounds = ubounds,
                        ubounds_rel = ubounds_rel,
                        dmu_eval = dmu_eval,
                        dmu_ref = dmu_ref
    )

    return(modelOutput)

  }
