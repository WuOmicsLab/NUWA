#' A portal function running xCell after NUWAms analysis.
#'
#' Run NUWAms and xCell algorithm with a marker list of interest.
#'
#' @param expr see the same argument of NUWA.cibersort.
#' @param marker_list a list, whose names are cellular populations' names
#' and elements are character vectors of markers (HUGO symbols), default is
#' xCell64.
#' @param ... additional arguments passed to the NUWAms() function
#'
#' @return see NUWA.cibersort.
#' @export
#'
#' @examples
#' expr <- cptacDatasets$brca[, 1:5]
#' res_nuwa <- NUWA.xcell(expr, marker_list = NULL)
#' res_nuwa <- NUWA.xcell(expr, marker_list = my_markers)
NUWA.xcell <- function(expr, marker_list = NULL, ...) {
    if(!require("xCell", quietly = TRUE)) {
        remotes::install_github('dviraran/xCell')
    } else {
        require(xCell)
    }
    if (!is.matrix(expr)) {
        stop("Parameter 'expr' should be a matrix")
    }

    if (is.null(marker_list)) {
        marker_list = xCell::xCell.data$signatures
        markers <- unlist(lapply(marker_list, function(x) x@geneIds))
    } else {
        markers <- unlist(marker_list)
    }


    args <- list(...)
    args <- modifyList(args, list(expr = expr, markers=unique(markers)), keep.null = T)

    res <- do.call(NUWAms, args)
    expr_impute <- res$finalExpr
    predVsTruth <- res$predVsTruth
    prop <- xCell::xCellAnalysis(expr_impute, signatures = marker_list,
                                 genes = rownames(expr_impute))
    return(list(proportion = t(prop), finalExpr = expr_impute, predVsTruth = predVsTruth))
}
