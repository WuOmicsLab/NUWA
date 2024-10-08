#' Immune cell types deconvolution using expression dataset.
#'
#' This function integrates deconvolution results of three algorithm-signature combinations selected from our benchmark analysis, and provides the relative proportions of six immune cell types in mixture samples. See the NUWA manuscript for more details.
#'
#' @param expr a numeric matrix of expression profiles for bulk tissue samples, with HUGO gene symbols as rownames and sample identifiers as colnames. Data must be non-logarithm scale.
#' @param cibersortPath a string specifying the path of CIBERSORT R script, CIBERSORT is only freely available for academic users, please register on https://cibersort.stanford.edu, and download the CIBERSORT source script.
#' 
#' @param BCIC_min_marker_num a positive integer, indicating the minimal number of BCIC markers needed to run EPIC. Default is 6.
#' @param LM6_min_marker_num a positive integer, indicating the minimal number of LM6 markers needed to run CIBERSORT. Default is 6.
#' @param LM22_min_marker_num a positive integer, indicating the minimal number of LM22 markers needed to run CIBERSORT. Default is 6.
#' 
#' @param quantileNorm_cibersort logical, indicating whether quantile normalization will be performed in CIBERSORT analysis. Only set FALSE for RNA-seq data as recommended on the CIBERSORT website. Default is TRUE.
#' @param protein logical, set TRUE for proteomic expression data. If TRUE, signature matrix including 118 markers (union of BCIC and TIC markers) will be used for EPIC analysis. If FALSE, the BCIC markers (n = 65) will be used. Default is TRUE.
#'
#' @return A list containing:\describe{
#'  \item{\code{proportion}}{ a matrix, the first column is the cell type name, and the remaining columns (one sample per column) are the relative proportion of mRNA or protein coming from the six immune cell types (B, CD4 T, CD8 T, monocyte/macrophage, NK and neutrophils cells).}
#'  \item{\code{mergedProp}}{a list containing three merged and sum-to-one normalized proportion matrices predicted by CIBERSORT-LM22, CIBERSORT-LM6 and EPIC-BCIC.}
#'  \item{\code{rawRes}}{a list containing the raw proportion matrices generated by CIBERSORT-LM22, CIBERSORT-LM6 and EPIC-BCIC.}
#'  \item{\code{usedComb}}{a string vector showing the combinations (from CIBERSORT-LM22, CIBERSORT-LM6 and EPIC-BCIC) used to generate the ensembled prediction of proportions. Disqualification might be caused by insufficient markers.}
#' }
#' @export
#'
#' @examples
#'  # You need to provide path to CIBERSORT.R
#'  # cibersortPath = "<PATHTO>/CIBERSORT.R"
#'  my.nuwadec = NUWAeDeconv(expr=my.nuwams$finalExpr, cibersortPath= cibersortPath)

NUWAeDeconv <- function(expr, cibersortPath, BCIC_min_marker_num = 6,
                       LM6_min_marker_num = 6, LM22_min_marker_num = 6,
                       quantileNorm_cibersort=T, protein=T) {
    if (!check_cibersort(cibersortPath)){
        stop("Invalid cibersort path. Please download script from cibersort website (https://cibersort.stanford.edu).")
    }
    filename=NULL
    mix <- expr
    RNAseq <- !quantileNorm_cibersort
    listmode <- F
    if (is.list(mix) & !is.data.frame(mix)) {
        mixls <- lapply(mix, function(x) {x[x < 0] <- 0; x})
        listmode <- T # Don't use this mode.
        if (length(mixls) == 0) stop("Error: listmode error")
    } else {
        if (!isSingleString(mix) && !is.matrix(mix)) {
            stop("'expr' needs to be given as a matrix or string path")
        }
        if (isSingleString(mix)) {
            mix <- read.delim(mix, row.names = 1, check.names = F)
            mix <- data.matrix(mix)
        }
        mix[mix < 0] <- 0
    }
    if (!is.null(filename) && length(filename) != 2) {
        stop("'filename' needs to be given as
            a string vector with 2 items or NULL")
    }
    combs <- c("epic_bcic", "cibersort_lm6", "cibersort_lm22")
    combid <- seq_along(combs)
    names(combid) <- names(combs) <- combs
    if (protein) {
        RNAseq <- F
        ## 118 genes
        siggs <- union(EPIC::BRef$sigGenes, EPIC::TRef$sigGenes)
    } else {
        ## 65 genes
        siggs <- EPIC::BRef$sigGenes
    }
    ref <- EPIC::BRef
    ref$sigGenes <- siggs

    if (listmode) {
        use <- sapply(combs, function(x) x %in% names(mixls))
    } else {
        siggsls <- list(siggs, rownames(LM6), rownames(LM22))
        names(siggsls) <- combs
        gs_mix <- rownames(mix)
        mar_num <- c(BCIC_min_marker_num,
                     LM6_min_marker_num, LM22_min_marker_num)
        use <- sapply(combid, function(i) {
            gs <- intersect(gs_mix, siggsls[[i]])
            return(length(gs) >= mar_num[i])
        })
    }
    if (listmode) mix <- mixls[["epic_bcic"]]
    epicRes0 <- tryCatch({
        options(warn = -1)
        epicgs <- intersect(rownames(EPIC::BRef$refProfiles), rownames(mix))
        if (length(epicgs) < 2e3) scale_expr <- F else scale_expr <- T
        a <- EPIC::EPIC(mix, ref, scaleExprs = scale_expr)
        options(warn = 0)
        a
    }, error = function(e) NULL)
    excluEpic <- identical(epicRes0, NULL)
    epicRes <- if (excluEpic) NULL else epicRes0[["mRNAProportions"]]
    qn <- !RNAseq
    if (listmode) mix <- mixls[["cibersort_lm6"]]
    cib_lm6_res <- tryCatch(CIBERSORT(LM6, mix, QN = qn),
                            error = function(e) NULL)
    if (listmode) mix <- mixls[["cibersort_lm22"]]
    cib_lm22_res <- tryCatch(CIBERSORT(LM22, mix, QN = qn),
                             error = function(e) NULL)
    resls <- list(epicRes, cib_lm6_res, cib_lm22_res)
    isNull <- sapply(resls, function(x) identical(x, NULL))
    use <- use & (!isNull)
    names(resls) <- combs
    if (sum(use) == 0) {
        cat("The number of genes included in both mixture and signatures:\n")
        marknum <- sapply(siggsls, function(x) length(intersect(gs_mix, x)))
        print(marknum)
        stop("There is no combination qualified. You can loose the
        restrictions on marker gene number based on the above infomation.")
    }
    resls_good <- lapply(resls[!isNull], function(x) x)
    celllist <- list(
        B = c("B cells naive", "B cells memory", "Bcells", "B cells"),
        CD4 = c("CD4_Tcells", "CD4 T cell", "CD4 T cells", "CD4.T.cells",
                "T cells CD4 naive", "T cells CD4 memory resting",
                "T cells CD4 memory activated",
                "T cells regulatory (Tregs)", "T cells follicular helper"),
        CD8 = c("T cells CD8", "CD8_Tcells", "CD8 T cells",
                "CD8 T cell", "CD8.T.cells"),
        NK = c("NKcells", "NK cells", "NK cell",
               "NK.cells", "NK cells activated"),
        Mono_Macro = c("Monocytes", "Macrophages M0",
                       "Macrophages M1", "Macrophages M2"),
        Neutro = c("Neutrophil", "Neutrophils")
    )
    resls_merged <- lapply(resls_good, mergecell, celllist = celllist)
    celllist1 <- celllist["CD4" != names(celllist)]
    resls_use <- resls_good[use[names(resls_good)]]
    resls5c <- lapply(resls_use, mergecell, celllist = celllist1)
    res_array5c <- abind::abind(resls5c, along = 3)
    res5c <- apply(res_array5c, c(1, 2), mean, na.rm = T)
    res <- resls_merged[["cibersort_lm22"]]
    res[, -2] <- (1 - res[, 2]) * res5c
    mRNAProp <- res
    if (!protein) {
        renormalize_vector <- c(
            "B" = 0.4016,
            "CD4" = 0.3952,
            "CD8" = 0.3952,
            "NK" = 0.4396,
            "Mono_Macro" = 1.4196,
            "Neutro" = 0.1300)
        res <- t(t(res) / renormalize_vector)
        res <- res / rowSums(res)
        cellProp <- res
    } else cellProp <- NULL
    fres <- list(
        proportion = mRNAProp,
        # cellProp = cellProp,
        mergedProp = resls_merged,
        rawRes = resls,
        usedComb = use
    )

    if (!is.null(filename)) {
        for (i in 1:2) {
            res <- fres[[i]]
            fname <- filename[i]
            if (!is.na(fname)) {
                res1 <- data.frame(SampleId = rownames(res),
                                   res, check.names = F)
                write.table(res1, fname, row.names = F, sep = "\t", quote = F)
            }
        }
    }
    return(fres)
}
