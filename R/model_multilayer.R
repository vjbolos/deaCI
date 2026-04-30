#' @title Multi-layer DEA-CI model
#'
#' @description Multi-layer DEA model for composite indicators, according to Shen at al. (2013).
#'
#' @usage model_multilayer(data_indicators,
#'                  layer_list = NULL,
#'                  hierarchy_tree = NULL,
#'                  sharebounds = NULL,
#'                  ubounds = NULL,
#'                  ubounds_rel = NULL,
#'                  wrange = NULL,
#'                  dmu_eval = NULL,
#'                  dmu_ref = NULL,
#'                  returnlp = FALSE)
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
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model and parameters \code{data_indicators}, \code{layer_list},
#'  \code{sharebounds}, \code{ubounds}, \code{ubounds_rel}, \code{wrange}, \code{dmu_eval}
#'  and \code{dmu_ref}.
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
#' # Example 1
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
#' result <- model_multilayer(data_indicators = Shen2013,
#'                            layer_list = layer_list)
#' CI(result)
#'
#' # Example 2
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
#' result <- model_multilayer(data_indicators = Shen2013,
#'                            hierarchy_tree = hierarchy_tree)
#' CI(result)
#'
#' # Example 3
#' # Weight bounds are personalized in layer 2.
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
#' wrange <- list(
#'   layer1 = NULL,
#'   layer2 = rbind(rep(0.4, 6),
#'                  rep(0.6, 6))
#' )
#' result <- model_multilayer(data_indicators = Shen2013,
#'                            layer_list = layer_list,
#'                            wrange = wrange)
#' CI(result)
#'
#'
#' @seealso \code{\link{model_multilayer_super}}, \code{\link{CI}}, \code{\link{scores_multilayer}},
#' \code{\link{model_basicCI}}
#'
#' @import lpSolve data.tree
#'
#' @export

model_multilayer <-
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

    if (!is.null(layer_list)) {
      # Reordenar para que la primera layer sea la que más items tiene
      longitudes <- lengths(layer_list)
      hierarchy <- layer_list[order(-longitudes)]
      for (i in 1:length(hierarchy)) {
        idxhijos <- 1
        if (i == 1) {
          nombres <- Inames
        } else {
          nombres <- names(hierarchy[[i - 1]])
        }
        for (j in 1:length(hierarchy[[i]])) {
          names(hierarchy[[i]][[j]]) <- nombres[idxhijos:(idxhijos + length(hierarchy[[i]][[j]]) - 1)]
          idxhijos <- idxhijos + length(hierarchy[[i]][[j]])
        }
      }

    } else if (!is.null(hierarchy_tree)) {
      hierarchy <- Node_to_layerlist(hierarchy_tree)
    } else {
      stop("Parameters layer_list or hierarchy_tree needed.")
    }

    nlayers <- length(hierarchy) + 1     # Número de capas
    s <- rep(0, nlayers)
    s[1] <- nI  # Número de indicadores

    for (k in 1:(nlayers - 1)) {
      if (any(sort(unlist(hierarchy[[k]])) != 1:s[k])) {
        stop("Hierarchy tree list not valid.")
      }
      s[k + 1] <- length(hierarchy[[k]]) # Número de categorías en cada capa
    }

    # --- A[[k]][[f]] Indicadores dentro de la categoría f de la capa k ---
    A <- vector("list", nlayers)
    names(A) <- paste0("layer", seq_len(nlayers))
    A[[1]] <- lapply(seq_len(nI), function(i) i)
    names(A[[1]]) <- Inames
    for (k in 2:nlayers) {
      A[[k]] <- hierarchy[[k - 1]]
      for (f in 1:s[k]) {
        A[[k]][[f]] <- unlist(A[[k - 1]][A[[k]][[f]]])
        names(A[[k]][[f]]) <- Inames[A[[k]][[f]]]
      }
    }

    # --- RESTRICCIONES DE PESOS ---

    if (!is.null(ubounds)) {

      sharebounds <- NULL
      ubounds_rel <- NULL
      if (length(ubounds) == 2) {
        ubounds <- matrix(c(ubounds[1], ubounds[2]), nrow = 2, ncol = s[nlayers])
      } else if ((nrow(ubounds) != 2) || (ncol(ubounds) != s[nlayers])) {
        stop("Parameter ubounds incorrect.")
      }
      rownames(ubounds) <- c("min", "max")
      colnames(ubounds) <- names(A[[nlayers]])

    } else if (!is.null(ubounds_rel)) {

      sharebounds <- NULL
      if (length(ubounds_rel) == 2) {
        ubounds_rel <- matrix(c(ubounds_rel[1], ubounds_rel[2]), nrow = 2, ncol = s[nlayers])
      } else if ((nrow(ubounds_rel) != 2) || (ncol(ubounds_rel) != s[nlayers])) {
        stop("Parameter ubounds_rel incorrect.")
      }
      rownames(ubounds_rel) <- c("min", "max")
      colnames(ubounds_rel) <- names(A[[nlayers]])

    } else {

      if (is.null(sharebounds)) {
        if (s[nlayers] > 2) {
          shareboundsmax <- 1 / (s[nlayers] - 1)
          shareboundsmin <- shareboundsmax / 5
        } else {
          shareboundsmin <- 0.2
          shareboundsmax <- 0.8
        }
        sharebounds <- matrix(c(shareboundsmin, shareboundsmax), nrow = 2, ncol = s[nlayers])
      } else if (length(sharebounds) == 2) {
        sharebounds <- matrix(c(sharebounds[1], sharebounds[2]), nrow = 2, ncol = s[nlayers])
      } else if ((nrow(sharebounds) != 2) || (ncol(sharebounds) != s[nlayers])) {
        stop("Parameter sharebounds incorrect.")
      }
      rownames(sharebounds) <- c("min", "max")
      colnames(sharebounds) <- names(A[[nlayers]])

    }

    wrange_general <- c(0.8, 1.2)
    if (is.null(wrange)) {
      wrange <- rep(list(NULL), nlayers - 1)
    } else if (!is.list(wrange) && (length(wrange) == 2)) {
      if (wrange[1] < 0 || wrange[1] > 1 || wrange[2] < 1) {
        stop("Parameter wrange incorrect: wrange[1] must be between 0 and 1, and wrange[2] must be greater than 1.")
      }
      wrange_general <- wrange
      wrange <- rep(list(NULL), nlayers - 1)
    } else if (is.list(wrange) && (length(wrange) != nlayers - 1)) {
      stop("Parameter wrange incorrect. It must be a vector of length 2 or a list of length (number of layers - 1).")
    }

    # Generate wrange if an element is NULL
    for (i in 1:(nlayers - 1)) {
      if (is.null(wrange[[i]])) {
        wrange[[i]] <- matrix(0, nrow = 2, ncol = s[i])
        rownames(wrange[[i]]) <- c("min", "max")
        colnames(wrange[[i]]) <- names(A[[i]])
        aux_ini <- 1
        for (f in 1:length(hierarchy[[i]])) {
          nif <- length(hierarchy[[i]][[f]])
          aux_fin <- aux_ini + nif - 1
          wrange[[i]][1, aux_ini:aux_fin] <- rep(wrange_general[1] / nif, nif)
          wrange[[i]][2, aux_ini:aux_fin] <- rep(wrange_general[2] / nif, nif)
          aux_ini <- aux_fin + 1
        }
      }
    }
    names(wrange) <- names(A)[-nlayers]

    # Generate H

    H <- hierarchy
    for (i in 1:(nlayers - 1)) {
      aux_ini <- 1
      for (f in 1:length(H[[i]])) {
        nif <- length(hierarchy[[i]][[f]])
        aux_fin <- aux_ini + nif - 1
        H[[i]][[f]] <- list(
          min = wrange[[i]][1, aux_ini:aux_fin],
          max = wrange[[i]][2, aux_ini:aux_fin]
        )
        aux_ini <- aux_fin + 1
      }
    }

    # --- Modelo Multilayer ---

    obj <- "max"

    # Constraints of 1st block
    f.con.1 <- Y[dmu_ref, ]
    f.dir.1 <- rep("<=", ndr)
    f.rhs.1 <- rep(1, ndr)

    # Constraint of 2nd block (rangos de pesos)

    f.con.min <- NULL
    f.con.max <- NULL
    aux <- rep(0, nI)
    for (i in 1:(nlayers - 1)) {
      k <- i + 1
      fw <- 1
      for (f in 1:s[k]) {
        nif <- length(hierarchy[[i]][[f]])
        if (nif > 1) {
          for (j in 1:nif) {
            aux1 <- aux
            aux1[A[[i]][[fw]]] <- 1
            aux2 <- aux
            aux2[A[[k]][[f]]] <- H[[i]][[f]]$min[j]
            f.con.min <- rbind(f.con.min, aux1 - aux2)
            aux2 <- aux
            aux2[A[[k]][[f]]] <- H[[i]][[f]]$max[j]
            f.con.max <- rbind(f.con.max, aux1 - aux2)
            fw <- fw + 1
          }
        } else {
          fw <- fw + 1
        }
      }
    }
    ncon <- nrow(f.con.max)
    f.dir.min <- rep(">=", ncon)
    f.dir.max <- rep("<=", ncon)
    f.rhs.min <- rep(0, ncon)
    f.rhs.max <- f.rhs.min

    DMU <- vector("list", nde)
    names(DMU) <- dmunames[dmu_eval]

    for (dmu_idx in 1:nde) {

      ii <- dmu_eval[dmu_idx]

      f.obj <- Y[ii, ]   # Función objetivo

      # Constraints of 3rd Block (Shares: cotas para u_hat*y_0)

      f.con.umin <- NULL
      f.con.umax <- NULL
      f.dir.umin <- rep(">=", s[nlayers])
      f.dir.umax <- rep("<=", s[nlayers])

      if (!is.null(ubounds)) {

        f.rhs.umin <- NULL
        f.rhs.umax <- NULL
        for (f in 1:s[nlayers]) {
          aux1 <- aux
          aux1[A[[nlayers]][[f]]] <- 1
          f.con.umin <- rbind(f.con.umin, aux1)
          f.rhs.umin <- c(f.rhs.umin, ubounds[1, f])
          f.rhs.umax <- c(f.rhs.umax, ubounds[2, f])
        }
        f.con.umax <- f.con.umin

      } else if (!is.null(ubounds_rel)) {

        for (f in 1:s[nlayers]) {
          aux21 <- rep(ubounds_rel[1, f], nI)
          aux22 <- rep(ubounds_rel[2, f], nI)
          aux1 <- aux
          aux1[A[[nlayers]][[f]]] <- 1
          f.con.umin <- rbind(f.con.umin, aux1 - aux21)
          f.con.umax <- rbind(f.con.umax, aux1 - aux22)
        }
        f.rhs.umin <- rep(0, s[nlayers])
        f.rhs.umax <- f.rhs.umin

      } else {

        for (f in 1:s[nlayers]) {
          aux21 <- rep(sharebounds[1, f], nI) * Y[ii, ]
          aux22 <- rep(sharebounds[2, f], nI) * Y[ii, ]
          aux1 <- aux
          aux1[A[[nlayers]][[f]]] <- 1
          aux1 <- aux1 * Y[ii, ]
          f.con.umin <- rbind(f.con.umin, aux1 - aux21)
          f.con.umax <- rbind(f.con.umax, aux1 - aux22)
        }
        f.rhs.umin <- rep(0, s[nlayers])
        f.rhs.umax <- f.rhs.umin

      }

      # All Constraints

      f.con <- rbind(f.con.1, f.con.min, f.con.max, f.con.umin, f.con.umax)
      f.dir <- c(f.dir.1, f.dir.min, f.dir.max, f.dir.umin, f.dir.umax)
      f.rhs <- c(f.rhs.1, f.rhs.min, f.rhs.max, f.rhs.umin, f.rhs.umax)

      if (returnlp) {
        u_hat <- rep(0,nI)
        names(u_hat) <- Inames
        var = list(u_hat = u_hat)
        DMU[[dmu_idx]] <- list(direction = obj, objective.in = f.obj, const.mat = f.con,
                               const.dir = f.dir, const.rhs = f.rhs, var = var)
      } else {

        res <- lp(obj, f.obj, f.con, f.dir, f.rhs)

        if (res$status == 0) {
          CI <- res$objval
          u_hat <- res$solution
          names(u_hat) <- Inames
          u <- rep(0, s[nlayers])
          for (f in 1:s[nlayers]) {
            u[f] <- sum(u_hat[A[[nlayers]][[f]]])
          }
          names(u) <- names(hierarchy[[nlayers - 1]])
          weights <- A[1:(nlayers - 1)]
          for (i in 1:(nlayers - 1)) {
            k <- i + 1
            fw <- 1
            weights[[i]] <- rep(0, s[i])
            for (f in 1:s[k]) {
              nif <- length(hierarchy[[i]][[f]])
              if (nif > 1) {
                for (j in 1:nif) {
                  weights[[i]][fw] <- sum(u_hat[A[[i]][[fw]]]) / sum(u_hat[A[[k]][[f]]])
                  fw <- fw + 1
                }
              } else {
                weights[[i]][fw] <- 1
                fw <- fw + 1
              }
            }
          }
          names(weights[[1]]) <- Inames
          if (nlayers > 2) {
            for (i in 2:(nlayers - 1)) {
              names(weights[[i]]) <- names(hierarchy[[i - 1]])
            }
          }
        } else {
          CI <- NA
          u_hat <- NA
          u <- NA
          weights <- NA
        }
        DMU[[dmu_idx]] <- list(CI = CI, u_hat = u_hat, u = u, weights = weights)
      }

    }

    modelOutput <- list(modelname = "multilayer",
                        DMU = DMU,
                        data_indicators = data_indicators,
                        layer_list = hierarchy,
                        sharebounds = sharebounds,
                        ubounds = ubounds,
                        ubounds_rel = ubounds_rel,
                        wrange = wrange,
                        dmu_eval = dmu_eval,
                        dmu_ref = dmu_ref
    )

    return(modelOutput)

  }
