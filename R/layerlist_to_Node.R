#' @title Layer list to Node
#'
#' @description Transforms a list
#' with the categories of each layer (from layer 2 upwards) and their subcategories
#' (as it is passed as parameter \code{layer_list} in function \code{model_multilayer}),
#' into a \code{Node} structure from package \pkg{data.tree}
#' (representing the hierarchical tree of categories and indicators).
#'
#' @usage layerlist_to_Node(layer_list)
#'
#' @param layer_list A list with the categories of each layer (from layer 2 upwards) and their subcategories,
#' as it is passed as parameter \code{layer_list} in function \code{model_multilayer}.
#'
#' @returns A \code{Node} structure from package
#' \pkg{data.tree}, representing the hierarchical tree of categories and indicators.
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
#'
#' # Replication of results in Shen et al. (2013) using hierarchy_tree.
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
#'
#' hierarchy_tree <- layerlist_to_Node(layer_list)
#' print(hierarchy_tree)
#' result <- model_multilayer(data_indicators = Shen2013, hierarchy_tree = hierarchy_tree)
#' CI(result)
#'
#' @seealso \code{\link{Node_to_layerlist}}, \code{\link{model_multilayer}}, \code{\link{model_multilayer_super}}
#'
#' @import lpSolve data.tree stats
#'
#' @export

layerlist_to_Node <-
  function(layer_list) {

    # Number of indicators
    nl2 <- length(layer_list$layer2)
    ni <- max(layer_list$layer2[[nl2]])

    # add layer1 with indicators
    layer1 <- setNames(
      replicate(ni, list(), simplify = FALSE),
      paste0("I", 1:ni)
    )
    layer_list$layer1 <- layer1
    layer_list <- layer_list[order(names(layer_list))]

    # construcción de capas

    nlayers <- length(layer_list)

    for (i in 2:nlayers) {
      # añado layer(i-1) a layer(i)
      nl <- length(layer_list[[i]])
      for (j in 1:nl) {
        layer_list[[i]][[j]] <- layer_list[[i - 1]][layer_list[[i]][[j]]]
      }
    }

    nodelist <- layer_list[[nlayers]]

    result_tree <- as.Node(nodelist)
    result_tree$name <- "CI"
    return(result_tree)

  }
