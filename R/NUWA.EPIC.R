#' Built-in NUWA analysis
#'
#' run NUWAms and EPIC algorithm with interested signature matrix
#'
#' @param expr see the same argument in NUWA.cibersort.
#' @param signature_matrix see the same argument in NUWA.cibersort.
#'
#' @return see NUWA.cibersort.
#' @export
#'
#' @examples
#' res_nuwa <- NUWAms.EPIC(expr = raw_expr, signature_matrix = BCIC)
#' res_nuwa <- NUWAms.EPIC(expr = raw_expr, signature_matrix = TIC)
#' res_nuwa <- NUWAms.EPIC(expr = raw_expr, signature_matrix = my_signature_matrix)
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
    predVsTruth <- res$predVsTruth
    ref <- list(refProfiles = signature_matrix, sigGenes = rownames(signature_matrix))
    prop <- EPIC::EPIC(expr_impute, ref)
    return(list(proportion = prop, mixture_impute = expr_impute, predVsTruth = predVsTruth))
}
