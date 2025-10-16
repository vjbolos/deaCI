#' @title CI scores
#'
#' @description Extract the composite indicator scores (optimal objective values) of the evaluated
#' DMUs from a multi-layer or basic DEA-CI solution returned by a \code{model_*} function.
#'
#' @usage CI(x)
#'
#' @param x List returned by a \code{model_*} function.
#'
#' @returns An array with CI scores of evaluated DMUs.
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
#'
#' @seealso \code{\link{model_multilayer}}, \code{\link{model_basicCI}}, \code{\link{model_multilayer_super}},
#'  \code{\link{model_basicCI_super}}, \code{\link{scores_multilayer}}
#'
#' @export

CI <-
  function(x) {

    result <- x

    if (result$modelname %in% c("multilayer", "basicCI", "multilayer_super", "basicCI_super")) {
      resCI <- unlist(lapply(result$DMU, function(x)
        x$CI))
    } else {
      stop("A list returned by function model_multilayer or model_basicCI is needed.")
    }

    return(resCI)

  }
