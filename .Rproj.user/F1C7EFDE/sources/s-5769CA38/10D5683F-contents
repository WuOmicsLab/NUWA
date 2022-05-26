# Data preprocessing
preprocess <- function(expr, batchInfo = NULL, thre = 50, nsams = 10, use.sams = F,
                       idcorr = F) {
    if (any(duplicated(colnames(expr)))) {
        stop("Sample ids should be unique.")
    }
    if (any(duplicated(rownames(expr)))) {
        stop("Gene symbols should be unique.")
    }

    # log2 trans
    if (!any(expr < 0, na.rm = T)) {
        expr <- log(expr, 2)
    }
    # geneid correct

    if (idcorr) {
        cat("Gene id correction ===============\n")
        print(dim(expr))
        expr <- geneidcorrect(expr)
        print(dim(expr))
        cat("Gene id correction finished ===============\n")
    }

    ## quantile normalization
    quantile_by_batch <- T
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
    if (quantile_by_batch & have_batch_info) {
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
    # 2 expr
    expr <- 2^expr.qn
    expr <- data.matrix(expr)
    # completeness >=50%
    expr <- completeness.ccb2(expr, threshold = thre, nsams = nsams, use.sams = use.sams)
    # cat("completeness 50%: remove ", nrow(expr.qn) - nrow(expr), " proteins\n")
    # remove outliers:
    expr <- rm.outlier(expr, usecap = T)
    # # matrix-level rescale the value to 0~1
    # expr <- global.scale(expr)
    return(expr)
}

