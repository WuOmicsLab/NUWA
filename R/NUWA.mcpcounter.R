#' Built-in NUWA analysis
#'
#' run NUWAms and MCPcounter algorithm with interested marker list.
#'
#' @param expr see the same argument in NUWA.cibersort.
#' @param marker_list see the same argument in NUWA.xcell, default is MCPcounter markers.
#'
#' @return see NUWA.cibersort.
#' @export
#'
#' @examples
#' res_nuwa <- NUWA.mcpcounter(expr = raw_expr, marker_list = NULL)
#' res_nuwa <- NUWA.mcpcounter(expr = raw_expr, marker_list = my_markers)
NUWA.mcpcounter <- function(expr, marker_list = NULL) {

    if (!is.matrix(expr)) {
        stop("Parameter 'expr' should be a matrix")
    }

    if (is.null(marker_list)) marker_list = mcpcounter_marker_list


    nw <- buildNetwork(markers = unique(unlist(marker_list)))
    res <- NUWAms(expr, network = nw)
    expr_impute <- res$finalExpr
    predVsTruth <- res$predVsTruth
    prop <- MCPcounter::appendSignatures(expr_impute, markers = marker_list)
    return(list(proportion = prop, mixture_impute = expr_impute, predVsTruth = predVsTruth))
}

