#' @title Radial super-efficiency multi-layer DEA-CI model
#'
#' @description Radial super-efficiency multi-layer DEA model for composite indicators,
#' according to Shen at al. (2013). Radial super-efficiency in the sense of Andersen
#' and Petersen (1993).
#'
#' @usage model_multilayer_super(data_indicators,
#'                        layer_list = NULL,
#'                        hierarchy_tree = NULL,
#'                        sharebounds = NULL,
#'                        ubounds = NULL,
#'                        ubounds_rel = NULL,
#'                        wrange = NULL,
#'                        dmu_eval = NULL,
#'                        dmu_ref = NULL,
#'                        returnlp = FALSE)
#'
#' @param data_indicators A data frame with the values of indicators of each DMU by rows. The names of
#'  the DMUs are the names of the rows. Optionally, the first column could contain the names of the DMUs.
#' @param layer_list A list with the categories of each layer (from layer 2 upwards)
#'  and their subcategories. Layer 1 (the bottom layer) is not included because it is assumed to always
#'  contain all the indicators, according to Shen et al. (2013).
#' @param hierarchy_tree Alternative to \code{layer_list}. It is a \code{Node} structure from package
#'  \pkg{data.tree}, representing the hierarchical tree of categories and indicators.
#'  Indicators are the leafs of the tree, not necessarily in the bottom layer. It is automatically
#'  transformed into a \code{layer_list} parameter. If \code{layer_list} is not \code{NULL}, it is ignored.
#' @param sharebounds It can be a vector of length \code{2} or a matrix of \code{2} rows and a number of
#'  columns equal to the number of indicators in the top layer. In the first case, it contains
#'  the lower and upper bounds (in proportion to 1) of the shares for all indicators in the top layer.
#'  The shares are the proportions that each indicator contributes to the efficiency score,
#'  according to Shen et al. (2013). In the second case, it contains the lower and upper bounds
#'  of the shares for each indicator in the top layer (minimums in the first row and maximums
#'  in the second row).
#'  If \code{sharebounds} is \code{NULL} (default), it is constructed automatically.
#'  For example, if there are 3 categories in the top layer, it takes the value \code{c(0.1, 0.5)},
#'  as in Shen et al. (2013).
#' @param ubounds It can be a vector of length \code{2} or a matrix of \code{2} rows and a number of
#'  columns equal to the number of indicators in the top layer. In the first case, it contains
#'  the lower and upper bounds of the weights u for all indicators in the top layer.
#'  In the second case, it contains the lower and upper bounds
#'  of the weights u for each indicator in the top layer (minimums in the first row and maximums
#'  in the second row). Note that the weights u do not add up to 1.
#' @param ubounds_rel It can be a vector of length \code{2} or a matrix of \code{2} rows and a number of
#'  columns equal to the number of indicators in the top layer. In the first case, it contains
#'  the lower and upper bounds of the weights u (relativized in proportion to 1) for all indicators in the top layer.
#'  In the second case, it contains the lower and upper bounds of the weights u (relativized
#'  in proportion to 1) for each indicator in the top layer (minimums in the first row and maximums
#'  in the second row). Note that the relativized weights sum to 1.
#' @param wrange It can be a vector of length \code{2} or a list of length \code{K - 1},
#'  where \code{K} is the number of layers.
#'
#'  In the first case, it contains the relative weight ranges for all layers (except
#'  the top layer) and categories. Weights vary in a range from \code{wrange[1]}
#'  to \code{wrange[2]} multiplied by the average weight in each category. Since the
#'  weights within a category must add up to 1, \code{wrange[1]} must be between 0 and 1,
#'  and \code{wrange[2]} must be greater than 1.
#'
#'  In the second case, the elements of the list are matrices of \code{2} rows containing
#'  the absolute weight ranges of the indicators within each layer (minimums in the first row
#'  and maximums in the second row), from layer \code{1} to layer \code{K - 1}.
#'  The number of columns must be equal to the number of indicators in the corresponding layer.
#'  If a matrix of the list (corresponding to the weight ranges of an entire layer) is \code{NULL}, then
#'  the weight ranges of this layer are constructed assuming that the relative weight ranges
#'  are \code{c(0.8, 1.2)}, as in Shen et al. (2013).
#'
#'  If \code{wrange} is \code{NULL} (default), it takes  the value \code{c(0.8, 1.2)},
#'  i.e. all layers have these relative weight bounds, as in Shen et al. (2013).
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation reference set.
#'  If \code{NULL} (default), all DMUs are considered.
#' @param returnlp Logical. If it is \code{TRUE}, it returns the linear problems (objective
#'  function and constraints).
#'
#' @returns A list of the results for the evaluated DMUs (\code{DMU} component),
#'  along with any other necessary information to replicate the results, such as the name of the model and
#'  parameters \code{data_indicators}, \code{layer_list}, \code{sharebounds}, \code{ubounds},
#'  \code{ubounds_rel}, \code{wrange}, \code{dmu_eval} and \code{dmu_ref}.
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
           sharebounds = NULL,
           ubounds = NULL,
           ubounds_rel = NULL,
           wrange = NULL,
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

  model_modelname <- "model_multilayer"

  DMU <- vector(mode = "list", length = nde)
  names(DMU) <- dmunames[dmu_eval]

  for (i in 1:nde) {

    ii <- dmu_eval[i]

    deasol <- do.call(model_modelname,
                      list(data_indicators = data_indicators,
                           layer_list = layer_list,
                           hierarchy_tree = hierarchy_tree,
                           sharebounds = sharebounds,
                           ubounds = ubounds,
                           ubounds_rel = ubounds_rel,
                           wrange = wrange,
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
                      sharebounds = deasol$sharebounds,
                      ubounds = deasol$ubounds,
                      ubounds_rel = deasol$ubounds_rel,
                      wrange = deasol$wrange,
                      dmu_eval = dmu_eval,
                      dmu_ref = dmu_ref
  )

  return(modelOutput)

}
