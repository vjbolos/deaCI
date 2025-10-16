#' @title Radial super-efficiency multi-layer DEA-CI model
#'
#' @description Radial super-efficiency multi-layer DEA model for composite indicators,
#' according to Shen at al. (2013). Radial super-efficiency in the sense of Andersen
#' and Petersen (1993).
#'
#' @usage model_multilayer_super(data_indicators,
#'                        layer_list = NULL,
#'                        hierarchy_tree = NULL,
#'                        ubounds = NULL,
#'                        wrange = NULL,
#'                        wrangematrix = NULL,
#'                        dmu_eval = NULL,
#'                        dmu_ref = NULL,
#'                        returnlp = FALSE)
#'
#' @param data_indicators A data frame with the values of indicators of each DMU by rows. The first column
#'  contains the names of the DMUs.
#' @param layer_list A list with the categories of each layer (from layer 2 upwards)
#'  and their subcategories. Layer 1 (the bottom layer) is not included because it is assumed to always
#'  contain all the indicators, according to Shen et al. (2013).
#' @param hierarchy_tree Alternative to \code{layer_list}. It is a \code{Node} structure from package
#'  \pkg{data.tree}, representing the hierarchical tree of categories and indicators.
#'  Indicators are the leafs of the tree, not necessarily in the bottom layer. It is automatically
#'  transformed into a \code{layer_list} parameter.
#' @param ubounds A vector with lower and upper relative bounds (in proportion to 1) for \code{u_hat * y_0}
#'  in the top layer, according to Shen et al. (2013). By default, it is constructed automatically.
#'  For example, if there are 3 categories in the top layer, it takes the value \code{c(0.1, 0.5)},
#'  as in Shen et al. (2013).
#' @param wrange A vector with two components. Weights vary in a range from \code{wrange[1]} to
#'  \code{wrange[2]} of the average weight in each category. Since the weights within a category must
#'  add up to 1, \code{wrange[1]} must be between 0 and 1, and \code{wrange[2]} must be greater than 1.
#'  By default, it takes the value \code{c(0.8, 1.2)}, as in Shen et al. (2013).
#' @param wrangematrix A matrix of size \code{2 x (K - 1)} (where \code{K} is the number of layers)
#'  with weight ranges from layer \code{1} to layer \code{K - 1}. It is an alternative to
#'  \code{wrange}, in the case of different weight ranges for each layer.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective
#'  function and constraints).
#'
#' @returns A list of the results for the evaluated DMUs (\code{DMU} component),
#'  along with any other necessary information to replicate the results, such as the name of the model and
#'  parameters \code{data_indicators}, \code{layer_list}, \code{ubounds}, \code{wrangematrix},
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
#' Andersen, P.; Petersen, N.C. (1993). "A procedure for ranking efficient units in
#' data envelopment analysis", Management Science, 39, 1261-1264.
#'
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
#' result <- model_multilayer_super(data_indicators = Shen2013, layer_list = layer_list)
#' CI(result)
#'
#' @seealso \code{\link{model_multilayer}}, \code{\link{CI}}, \code{\link{scores_multilayer}},
#' \code{\link{model_basicCI_super}}
#'
#' @export

model_multilayer_super <-
  function(data_indicators,
           layer_list = NULL,
           hierarchy_tree = NULL,
           ubounds = NULL,
           wrange = NULL,
           wrangematrix = NULL,
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

  model_modelname <- "model_multilayer"

  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]

  for (i in 1:nde) {

    ii <- dmu_eval[i]

    deasol <- do.call(model_modelname,
                      list(data_indicators = data_indicators,
                           layer_list = layer_list,
                           hierarchy_tree = hierarchy_tree,
                           ubounds = ubounds,
                           wrange = wrange,
                           wrangematrix = wrangematrix,
                           dmu_eval = ii,
                           dmu_ref = dmu_ref[dmu_ref != ii],
                           returnlp = returnlp
                           )
                      )

    DMU[[i]] <- deasol$DMU[[1]]

  }

  modelOutput <- list(modelname = "multilayer_super",
                      DMU = DMU,
                      data_indicators = data_indicators,
                      layer_list = deasol$layer_list,
                      ubounds = deasol$ubounds,
                      wrangematrix = deasol$wrangematrix,
                      dmu_eval = dmu_eval,
                      dmu_ref = dmu_ref
  )

  return(modelOutput)

}
