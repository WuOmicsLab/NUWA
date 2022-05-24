#' A portal function running CIBERSORT after NUWAms analysis.
#'
#' Run NUWAms and MCPcounter algorithm with a marker list of interest.
#'
#' @param expr see the same argument in NUWA.cibersort.
#' @param marker_list see the same argument in NUWA.xcell, default is MCPcounter markers.
#' @param ... additional arguments passed to the NUWAms() function
#'
#' @return see NUWA.cibersort.
#' @export
#'
#' @examples
#' expr <- cptacDatasets$brca[, 1:5]
#' res_nuwa <- NUWA.mcpcounter(expr, marker_list = NULL)
#' res_nuwa <- NUWA.mcpcounter(expr, marker_list = my_markers)
NUWA.mcpcounter <- function(expr, marker_list = NULL, ...) {

    if (!is.matrix(expr)) {
        stop("Parameter 'expr' should be a matrix")
    }

    if (is.null(marker_list)) marker_list = mcpcounter_marker_list


    nw <- buildNetwork(markers = unique(unlist(marker_list)))
    args <- list(...)
    args <- modifyList(args, list(expr = expr, network=nw), keep.null = T)

    res <- do.call(NUWAms, args)
    expr_impute <- res$finalExpr
    predVsTruth <- res$predVsTruth
    prop <- MCPcounter::appendSignatures(expr_impute, markers = marker_list)
    return(list(proportion = prop, finalExpr = expr_impute, predVsTruth = predVsTruth))
}

