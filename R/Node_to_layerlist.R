#' @title Node to layer list
#'
#' @description Transforms a \code{Node} structure from package \pkg{data.tree}
#' (representing the hierarchical tree of categories and indicators) into a list
#' with the categories of each layer (from layer 2 upwards) and their subcategories,
#' as it is passed as parameter \code{layer_list} in function \code{model_multilayer}.
#'
#' @usage Node_to_layerlist(hierarchy_tree)
#'
#' @param hierarchy_tree A \code{Node} structure from package
#' \pkg{data.tree}, representing the hierarchical tree of categories and indicators.
#'  Indicators are the leafs of the tree, not necessarily in the bottom layer.
#'
#' @returns A list with the categories of each layer (from layer 2 upwards) and their subcategories,
#' as it is passed as parameter \code{layer_list} in function \code{model_multilayer}.
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
#' library(data.tree)
#' hierarchy_tree <- Node$new("CI")
#' # Layer 3
#' Alcohol <- hierarchy_tree$AddChild("Alcohol")
#' Speed <- hierarchy_tree$AddChild("Speed")
#' Protective_System <- hierarchy_tree$AddChild("Protective_System")
#' # Layer 2
#' I1 <- Alcohol$AddChild("I1")
#' I2 <- Alcohol$AddChild("I2")
#' Mean_Speed <- Speed$AddChild("Mean_Speed")
#' Speed_Limit <- Speed$AddChild("Speed_Limit")
#' Seat_Belt <- Protective_System$AddChild("Seat_Belt")
#' Child_Restraint <- Protective_System$AddChild("Child_Restraint")
#' # Layer 1
#' Mean_Speed$AddChild("I3")
#' Mean_Speed$AddChild("I4")
#' Mean_Speed$AddChild("I5")
#' Speed_Limit$AddChild("I6")
#' Speed_Limit$AddChild("I7")
#' Speed_Limit$AddChild("I8")
#' Seat_Belt$AddChild("I9")
#' Seat_Belt$AddChild("I10")
#' Child_Restraint$AddChild("I11")
#' print(hierarchy_tree)
#'
#' layer_list <- Node_to_layerlist(hierarchy_tree)
#' result <- model_multilayer(data_indicators = Shen2013, layer_list = layer_list)
#' CI(result)
#'
#' @seealso \code{\link{layerlist_to_Node}}, \code{\link{model_multilayer}}, \code{\link{model_multilayer_super}}
#'
#' @import lpSolve data.tree
#'
#' @export

Node_to_layerlist <-
  function(hierarchy_tree) {

    # Llevar Indicadores a la profundidad máxima
    profundidad_max <- max(hierarchy_tree$Get("level"))
    hojas <- hierarchy_tree$leaves
    for (hoja in hojas) {
      nivel_actual <- hoja$level
      nodo_actual <- hoja
      while (nivel_actual < profundidad_max) {
        nuevo_hijo <- nodo_actual$AddChild(nodo_actual$name)
        nodo_actual <- nuevo_hijo
        nivel_actual <- nivel_actual + 1
      }
    }
    # Construir la lista hierarchy a partir del árbol
    nodos_por_capa <- lapply(2:profundidad_max, function(nivel) {
      hierarchy_tree$Get("name", filterFun = function(node) node$level == nivel)
    })
    nl <- length(nodos_por_capa)
    nombres_h <- NULL
    for (i in 1:nl) {
      nombres_h <- c(nombres_h, paste0("layer", nl - i + 1))
    }
    names(nodos_por_capa) <- nombres_h

    hierarchy <- vector("list", nl - 1)
    for (i in 1:(nl - 1)) {
      nci <- length(nodos_por_capa[[i]])
      kk <- vector("list", nci)
      names(kk) <- nodos_por_capa[[i]]
      idxhijos <- 1
      for (j in 1:length(nodos_por_capa[[i]])) {
        nodo_madre <- data.tree::FindNode(hierarchy_tree, nodos_por_capa[[i]][j])
        kk[[j]] <- idxhijos:(idxhijos + length(names(nodo_madre$children)) - 1)
        names(kk[[j]]) <- names(nodo_madre$children)
        idxhijos <- idxhijos + length(names(nodo_madre$children))
      }
      hierarchy[[i]] <- kk
    }
    names(hierarchy) <- nombres_h[1:(length(nombres_h) - 1)]
    hierarchy <- rev(hierarchy)

    return(hierarchy)

  }
