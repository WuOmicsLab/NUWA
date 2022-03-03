
#' Built-in NUWA analysis
#'
#' run NUWAms and CIBERSORT algorithm with interested signature matrix.
#'
#' @param expr a numeric matrix or data frame of expression profiles for
#' bulk tissue samples, with HUGO gene symbols as rownames and sample
#' identifiers as colnames. It can also be a string specifing the file
#' path of an expression matrix. Data must be non-logarithm scale.
#' @param cibersortPath a string specifying the path of CIBERSORT R script,
#' CIBERSORT is only freely available for academic users, please register
#' on https://cibersort.stanford.edu, and download the CIBERSORT source script.
#' @param signature_matrix a signature matrix (such as lm22, lm6, BCIC, TIC).
#'
#' @return The results of each built-in NUWA analysis function, is a list
#' containing an expression matrix with missing markers inferred,
#' two matrices using for compute recall, and a matrix including cell
#' fractions estimated by the algorithm used.
#' @export
#'
#' @examples
#' res_nuwa <- NUWAms.cibersort(expr = raw_expr, cibersortPath= ciberR, signature_matrix = LM22)
#' res_nuwa <- NUWAms.cibersort(expr = raw_expr, cibersortPath= ciberR, signature_matrix = LM6)
#' res_nuwa <- NUWAms.cibersort(expr = raw_expr, cibersortPath= ciberR, signature_matrix = my_signature_matrix)
NUWA.cibersort <- function(expr, signature_matrix, cibersortPath) {

    impute <- T
    buildNetworkArgs <- list()
    nuwamsArgs <- list()
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
    if (impute) {

        buildNetworkArgs0 <- list(markers = rownames(signature_matrix))
        buildNetworkArgs0 <- modifyList(buildNetworkArgs0, buildNetworkArgs)
        nw <- do.call(buildNetwork, buildNetworkArgs0)

        nuwamsArgs0 <- list(expr = mix, network = nw)
        nuwamsArgs0 <- modifyList(nuwamsArgs0, nuwamsArgs)
        mix0 <- mix
        res <- do.call(NUWAms, nuwamsArgs0)
        mix <- res$finalExpr
        predVsTruth <- res$predVsTruth
        # save(signature_matrix, mix0, mix, nw, file = "/home/pub/project/deconvolution/tcga/cptac_test/nuwa_package_development/res_xiergo_new.rda")
    }
    prop <- CIBERSORT(signature_matrix, mix, QN = T)
    if (!impute) {
        mix <- NULL
        predVsTruth <- NULL
    }
    return(list(proportion = prop, mixture_impute = mix, predVsTruth = predVsTruth))
}
