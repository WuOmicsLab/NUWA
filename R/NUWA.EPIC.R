#' A portal function running EPIC after NUWAms analysis.
#'
#' Run NUWAms and EPIC algorithm with a signature matrix of interest.
#'
#' @param expr see the same argument in NUWA.cibersort.
#' @param signature_matrix see the same argument in NUWA.cibersort.
#'
#' @return see NUWA.cibersort.
#' @export
#'
#' @examples
#' expr <- cptacDatasets$brca[, 1:5]
#' res_nuwa <- NUWA.EPIC(expr, signature_matrix = BCIC)
#' res_nuwa <- NUWA.EPIC(expr, signature_matrix = TIC)
#' res_nuwa <- NUWA.EPIC(expr, signature_matrix = my_signature_matrix)
NUWA.EPIC <- function(expr, signature_matrix) {
    if (!is.matrix(expr)) {
        stop("Parameter 'expr' should be a matrix")
    }
    if (!is.matrix(signature_matrix)) {
        stop("Parameter 'signature_matrix' should be given as matrix.")
    }
    nw <- buildNetwork(markers = rownames(signature_matrix))
    res <- NUWAms(expr, network = nw)
    expr_impute <- res$finalExpr
    print(str(expr))
    print(str(expr_impute))
    print(str(res))
    predVsTruth <- res$predVsTruth
    ref <- list(refProfiles = signature_matrix, sigGenes = rownames(signature_matrix))
    prop <- EPIC::EPIC(expr_impute, ref)$mRNAProportions
    return(list(proportion = prop, finalExpr = expr_impute, predVsTruth = predVsTruth))
}
