#' A portal function running MCPcounter after NUWAms analysis.
#'
#' Run NUWAms and MCPcounter algorithm with a marker list of interest.
#'
#' @param expr see the same argument in NUWA.cibersort.
#' @param marker_list a list, which names are cellular populations' names and elements are character vectors of markers (HUGO symbols). Default is MCPcounter markers.
#' @param ... additional arguments passed to the NUWAms() function
#'
#' @return see NUWA.cibersort.
#' @export
#'
#' @examples
#' expr <- CPTAC.6datasets$brca[, 1:5]
#' res_nuwa <- NUWA.mcpcounter(expr, marker_list = NULL)
#' res_nuwa <- NUWA.mcpcounter(expr, marker_list = my_markers)

NUWA.mcpcounter <- function(expr, marker_list = NULL, ...) {

    if (!is.matrix(expr)) {
        stop("Parameter 'expr' should be a matrix")
    }

    if (is.null(marker_list)) marker_list = MCPCOUNTER_MARKER_LIST


    args <- list(...)
    args <- modifyList(args, list(expr = expr, markers=unique(unlist(marker_list))), keep.null = T)

    res <- do.call(NUWAms, args)
    expr_impute <- res$finalExpr
    predVsTruth <- res$predVsTruth
    prop <- MCPcounter::appendSignatures(expr_impute, markers = marker_list)
    return(list(proportion = prop, finalExpr = expr_impute, predVsTruth = predVsTruth))
}

