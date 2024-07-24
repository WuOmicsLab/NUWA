#' A portal function running CIBERSORT after NUWAms analysis.
#'
#' Run NUWAms (missing markers inference) and CIBERSORT algorithm (deconvolution) with a signature matrix of interest.
#'
#' @param expr a numeric matrix of expression profiles for bulk tissue samples, with HUGO gene symbols as rownames and sample identifiers as colnames. Data must be non-logarithm scale.
#' @param cibersortPath a string specifying the path of CIBERSORT R script, CIBERSORT is only freely available for academic users, please register on https://cibersort.stanford.edu, and download the CIBERSORT source script.
#' @param signature_matrix a signature matrix, which is a numeric expression matrix of markers in cell types of interest, with HUGO gene symbols as rownames and cell type identifiers as colnames. Such as LM22, LM6, BCIC, TIC or user provided signature matrix.
#' @param ... additional arguments passed to the NUWAms() function
#'
#' @return The results of each built-in NUWA analysis function, is a list containing an expression matrix with missing markers inferred, two matrices used for recall analysis, and a matrix including cell fractions estimated by the algorithm used.
#' @export
#'
#' @examples
#' expr <- cptacDatasets$brca[, 1:5]
#' res_nuwa <- NUWA.cibersort(expr, cibersortPath = cibersortPath, signature_matrix = LM22)
#' res_nuwa <- NUWA.cibersort(expr, cibersortPath = cibersortPath, signature_matrix = NUWAp26)
#' res_nuwa <- NUWA.cibersort(expr, cibersortPath = cibersortPath, 
#'                          signature_matrix = my_signature_matrix)
NUWA.cibersort <- function(expr, signature_matrix, cibersortPath, ...) {
    if (!check_cibersort(cibersortPath)){
        stop("Invalid cibersort path. Please download script from cibersort website (https://cibersort.stanford.edu).")
    }
    if (! is.matrix(expr)) {
        stop("'expr' should be a matrix")
    }
    mix <- expr

    if (!is.matrix(signature_matrix)) {
        stop("Parameter signature_matrix should be given as matrix.")
    }
    if (!file.exists(cibersortPath)) {
        stop("Could not find file ", cibersortPath)
    }

    nuwamsArgs0 <- modifyList(list(...), list(expr = mix, markers = rownames(signature_matrix)), keep.null = T)
    res <- do.call(NUWAms, nuwamsArgs0)
    mix <- res$finalExpr
    predVsTruth <- res$predVsTruth

    prop <- CIBERSORT(signature_matrix, mix, QN = T)
    if (!impute) {
        mix <- NULL
        predVsTruth <- NULL
    }
    return(list(proportion = prop, finalExpr = mix, predVsTruth = predVsTruth))
}
