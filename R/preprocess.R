##
#' Preprocess the numeric matrix of expression profiles.
#'
#' This function serves as an interface for preprocessing raw data, encompassing tasks such as gene ID correction, quantile normalization, filtering out samples with excessive missing values, and eliminating outliers within each sample.
#'
#' @param expr a numeric matrix of expression profiles for bulk tissue samples, with HUGO gene symbols as rownames and sample identifiers as colnames.
#' @param idcorr logical. If set to TRUE, correction will be applied to gene IDs using the gene names provided at https://www.genenames.org. Default FALSE.
#' @param quantile_normalization logical. If TURE, quantile normalization will be performed, default TRUE.
#' 
#' @param batchInfo A vector or factor of the same length as the number of samples is provided to indicate which samples belong to which batch. The order of this vector should be consistent with the order of the samples. If the order is not consistent, a sample ID-named vector or factor should be used. If this information is provided, quantile normalization is performed within each respective batch. Default NULL.
#' 
#' @param use.sams logical. If TRUE, genes will be filtered based on the absolute count of non-NA values; otherwise, genes will be filtered based on the relative proportion of non-NA values. Default FALSE.
#' @param thre A threshold ranging from 0 to 100. when use.sams is False, genes with a non-NA value proportion below `thre\%` will be filtered out. Default 0.
#' @param nsams a cutoff. When use.sams is TRUE, genes with a count of non NA lower than nsams will be filtered out. Default 10.
#' @param usecap If set to TRUE, lower outliers are replaced with the 5th percentile for each sample, and upper outliers are replaced with the 95th percentile; otherwise, they are replaced with NA. Here low outliers are below Q1 - 1.5 * IQR, and high outliers are above Q3 + 1.5 * IQR. Default TRUE.
#' 
#' @param quantification_method The quantification method used in MS-proteomic profiling, "TMT/iTRAQ ratio" or "Label-free intensity". Default is "TMT/iTRAQ ratio".
#' @param logbase The log base of the expression matrix if the quantification values have been log-transformed, one of "No transformation" or "Log2". Default is "No transformation".
#' @return A numeric matrix of expression profiles after preprocessing.
#'
#' @export
#'
#' @examples
#' expr <- CPTAC.6datasets$brca[, 1:5]
#' res <- preprocess(expr,
#'             quantification_method = "TMT/iTRAQ ratio",
#'             logbase = "No transformation")

preprocess <- function(expr, idcorr = F, quantile_normalization = T, 
                       batchInfo = NULL, use.sams = F,
                       thre = 0, nsams = 10, usecap = T,
                       quantification_method=c("TMT/iTRAQ ratio", "Label-free intensity")[1],
                       logbase = c("No transformation", "Log2")[1]
                       ) {
    
    require(matrixStats)
    if (any(duplicated(colnames(expr)))) {
        stop("Sample ids should be unique.")
    }
    if (any(duplicated(rownames(expr)))) {
        stop("Gene symbols should be unique.")
    }

    # non-log trans
    if (any(expr < 0, na.rm = T) | logbase == "Log2") {
        expr <- 2^expr
    }

    # TMT-like trans for label-free data (non-log scale)
    if (quantification_method == "Label-free intensity") {
        ## row-level (within gene) ratio trans: exp_row/median_row
        v0 <- rowMedians(expr, na.rm=T)
        mat_r <- expr/v0

        ## column-level (within sample) ratio trans: exp_column/median_column
        v1 <- colMedians(mat_r, na.rm=T)
        expr <- t(apply(mat_r, 1, function(x) {x/v1}))
    }

    # log2 trans
    expr <- log(expr, 2)

    # geneid correct
    if (idcorr) {
        cat("Gene id correction ===============\n")
        print(dim(expr))
        expr <- geneidcorrect(expr)
        print(dim(expr))
        cat("Gene id correction finished ===============\n")
    }

    ## quantile normalization
    have_batch_info <- !is.null(batchInfo)
    if (have_batch_info) {
        if (!is.vector(batchInfo) | !(is.factor(batchInfo))) {
            stop("batchInfo should be a vector or factor.")
        }
        if (length(batchInfo) != ncol(expr)) {
            stop("The length of batchInfo should be equal to the number of samples in expr.")
        }
        if (!is.null(names(batchInfo))) {
            if (length(intersect(names(batchInfo), colnames(expr))) != length(colnames(expr))) {
                stop("The names of batchInfo should be the sample ids in expr.")
            }
            batchInfo <- batchInfo[match(colnames(expr), names(batchInfo))]
        }
    }

    ##
    if (quantile_normalization & have_batch_info) {
        id_batch_df <- data.frame(
            SampleID = colnames(expr),
            batch = batchInfo
        )
        # print(head(id_batch_df))
        batchs <- as.character(id_batch_df$batch)
        # expr0 <- expr
        for (b in unique(batchs)) {
            # print(b)
            expr[, batchs == b] <- preprocessCore::normalize.quantiles(as.matrix(expr[, batchs == b]), copy = TRUE)
        }

        # remove genes with NAs in all samples of one batch.
        not_na_mat <- sapply(
            unique(batchs),
            function(b) {
                apply(!is.na(expr[, batchs == b]), 1, sum)
            }
        )
        min_not_na <- apply(not_na_mat, 1, min)
        print(table(min_not_na))
        expr <- expr[min_not_na != 0, ]
        # remove batch effect
        modcombat <- model.matrix(~1, data = id_batch_df)
        expr <- sva::ComBat(dat = expr, batch = batchs, mod = modcombat, par.prior = TRUE, prior.plot = FALSE, mean.only = TRUE)
        expr.qn <- expr
    } else {
        expr.qn <- preprocessCore::normalize.quantiles(as.matrix(expr), copy = TRUE)
        rownames(expr.qn) <- rownames(expr)
        colnames(expr.qn) <- colnames(expr)
    }

    # debug(expr.qn,5); dim(expr.qn)
    # 2^expr
    expr <- 2^expr.qn
    expr <- data.matrix(expr)
    # completeness >= xx%
    expr <- completeness.ccb2(expr, threshold = thre, nsams = nsams, use.sams = use.sams)
    # cat("completeness xx%: remove ", nrow(expr.qn) - nrow(expr), " proteins\n")
    # remove outliers:
    expr <- rm.outlier(expr, usecap = T)
    # # matrix-level rescale the value to 0~1
    # expr <- global.scale(expr)
    return(expr)
}

