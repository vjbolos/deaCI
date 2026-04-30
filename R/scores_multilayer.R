#' @title Category scores in a multi-layer DEA-CI model
#'
#' @description Extract category scores of the evaluated DMUs from a multi-layer
#' DEA-CI solution returned by function \code{model_multilayer} or \code{model_multilayer_super}.
#'
#' @usage scores_multilayer(x)
#'
#' @param x List returned by function \code{model_multilayer} or or \code{model_multilayer_super}.
#'
#' @returns A list with the category scores of the evaluated DMUs.
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
#' layer_list <- list(
#'   layer2 = list(
#'     I1 = c(1L),
#'     I2 = c(2L),
#'     Mean_Speed = c(3L, 4L, 5L),
#'     Speed_Limit = c(6L, 7L, 8L),
#'     Seat_Belt = c(9L, 10L),
#'     Child_Restraint = c(11L)
#'   ),
#'   layer3 = list(
#'     Alcohol = c(1L, 2L),
#'     Speed = c(3L, 4L),
#'     Protective_System = c(5L, 6L)
#'   )
#' )
#' result <- model_multilayer(data_indicators = Shen2013, layer_list = layer_list)
#' CI(result)
#' category_scores <- scores_multilayer(result)
#'
#' @seealso \code{\link{model_multilayer}}, \code{\link{CI}}
#'
#' @export

scores_multilayer <-
  function(x) {

    if (!(x$modelname %in% c("multilayer", "multilayer_super"))) {
      stop("A list returned by function model_multilayer or model_multilayer_super is needed.")
    }

    dmu_eval <- x$dmu_eval
    Y_eval <- as.matrix(x$data_indicators[dmu_eval, ])

    nde <- length(x$DMU)
    scores <- vector(mode = "list", length = nde)
    names(scores) <- names(x$DMU)

    layer_list <- x$layer_list
    nll <- length(layer_list)
    for (j in 1:nde) {

      scores_dmu <- vector(mode = "list", length = nll)
      names(scores_dmu) <- names(layer_list)

      aux <- Y_eval[j, ]
      for (k in 1:nll) {
        wY <- x$DMU[[j]]$weights[[k]] * aux
        nlli <- length(layer_list[[k]])
        scores_dmu[[k]] <- rep(0, nlli)
        names(scores_dmu[[k]]) <- names(layer_list[[k]])
        for (i in 1: nlli) {
          scores_dmu[[k]][i] <- sum(wY[unlist(layer_list[[k]][i])])
        }
        aux <- scores_dmu[[k]]
      }

      scores[[j]] <- scores_dmu

    }

    return(scores)

  }


