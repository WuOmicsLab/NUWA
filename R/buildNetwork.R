
#' Build co-expression network for individual marker based on provided training datasets.
#'
#' This function builds individual co-expression network for each provided marker,
#' to select proteins with correlated expression relationship of a marker
#' (using samples with marker quantification) in the given training datasets
#' (two or more proteome datasets needed).
#'
#' @param trainsets a list containing the underlying datasets used for building the co-expression networks and for training regression models in the subsequent analysis. Each dataset in the list should be a numeric expression matrix (non-logarithm scale), with HUGO gene symbols as rownames and sample identifiers as colnames. When default (NULL) is applied, a list including CPTAC proteomic datasets of six cancer types (breast cancer, clear cell renal cell carcinoma, colon adenocarcinoma, endometrial carcinoma, gastric cancer, lung adenocarcinoma) will be used (CPTAC.6datasets). Users can provide customized datasets with more specific context, e.g., multiple datasets for a specific cancer type.
#'
#' @param markers a character vector of interesting proteins (termed as marker) that need to build co-expression networks. If NULL (default) provided, a pre-defined list MARKERS_Immune.1114 (1114 immune marker proteins collected from public studies) will be used. Another two pre-defined lists, including MARKER_FDA_Drug_Targets.812 (FDA-approved drug targets from HPA and Drugbank) and MARKERS_Immune.NUWAp26 (631 immune cell markers from proteomic signature matrix NUWAp26 we developed for 26 immune cell types) were also provided for convenience.
#'
#' @param preprocess logical. If TURE, each matrix in "trainsets" is preprocessed before network building. Default is TRUE. See the Methods section of the NUWA manuscript for more details.
#'
#' @param batchInfoList a list containing the batch information corresponding to the training datasets. The batchInfoList should be named using the names of trainsets, or the length of batchInfoList should be equal to the number of trainsets. Each element in the list gives the batch information of one trainset, which is a vector named by the sample identifiers of that trainset.
#'
#' @param nTr a positive integer, we only build network for markers existing in greater than or equal to "nTr" training datasets. Default is 2.
#'
#' @param corCutoff a positive numeric, specifying the absolute Pearson correlation coefficient threshold above which a co-expression will be declared between individual marker and other quantified protein. For each marker, we identify its "coherently correlated proteins", i.e., proteins with the same correlation coefficient sign in all training datasets, and with significant correlation (P < 0.05, absolute value of Pearson correlation coefficient greater than "corCutoff") in at least two training datasets. Default is 0.3.
#'
#' @param ncores a positive integer, indicating the number of cores used by this function. If the operating system is windows, then only one core will be used.
#'
#' @return A list containing:\describe{
#'  \item{\code{corr}}{a numeric data frame of Pearson correlation coefficients between markers and other proteins.}
#'  \item{\code{markers}}{a character vector, containing the markers with co-expression networks.}
#'  \item{\code{trainScaled}}{a list containing the preprocessed training datasets.}
#'  }
#' @export
#'
#' @examples
#' my.net <- buildNetwork(CPTAC.6datasets[1:5])
#' str(my.net)
#' 
buildNetwork <- function(trainsets = NULL, markers = NULL, preprocess = T,
                          batchInfoList = NULL, nTr = 2,
                            corCutoff = 0.3, ncores = 16) {

    if (is.null(markers)) markers <- MARKERS_Immune.1114
    input_markers <- markers
    if (is.null(trainsets)) {
        cat("Using six cptac datasets\n")
        if (nTr == 2 & corCutoff == 0.3 & (!preprocess)) {
            cat("Using pre-computed network\n")
            nw <- pruneNetwork(markers)
            # nw <- NETWORK
            # gps <- nw$genepair
            # gpls <- strsplit(gps, "~")
            # markers <- intersect(markers, NETWORK$markers)
            # gps <- gps[sapply(gpls, function(x) {any(x %in% markers)})]
            # nw$markers <- markers
            # nw$genepair <- gps
            # trs_scaled <- lapply(CPTAC.6datasets, function(exp) exp <- global.scale(exp)$y)
            # nw$trainScaled <- trs_scaled
            return(nw)
        } else {
            trainsets <- CPTAC.6datasets
        }
    }

    if (length(trainsets) < 3) {
        stop("There should be at least three training datasets in trainsets.")
    }
    if (is.null(names(trainsets))) {
        names(trainsets) <- paste0("Trainset_", 1:length(trainsets))
    }
    if (!is.null(batchInfoList)) {
        if (is.null(names(batchInfoList))) {
            if (length(batchInfoList) != length(trainsets)) {
                stop("The names of batchInfoList should be specified, or the length of batchInfoList should be equal to the number of trainsets")
            } else {
                names(batchInfoList) <- names(trainsets)
            }
        }
        batchInfoList <- batchInfoList[names(trainsets)]
    }
    ds_train <- names(trainsets)
    # caculate correlation between stable-marker and genes -------------------------------------------
    overlap.genes <- list()
    EXP <- list()
    cat("Dataset preprocessing -----------------------\n")
    for (I in 1:length(trainsets)) {
        cat("Dataset", I, "/", length(trainsets))
        exp <- trainsets[[I]]
        if (preprocess) {
            batchInfo0 <- batchInfoList[[I]]
            exp <- preprocess(expr = exp, batchInfo = batchInfo0)
        }
        # matrix-level rescale the value to 0~1
        exp <- global.scale(exp)$y
        cat("------- Finished\n")
        EXP[[I]] <- exp
        names(EXP)[I] <- ds_train[I]
        overlap.genes[[I]] <- intersect(rownames(exp), markers)
    }
    cat("\n\n")
    overlap <- unlist(overlap.genes)
    overlap <- as.data.frame(table(overlap))
    overlap <- as.character(overlap[overlap[, 2] >= nTr, 1])
    cordfls <- list()
    cat("Network construction -----------------------\n")
    for (I in 1:length(EXP)) {
        cat("Dataset", I, "/", length(EXP))
        exp <- EXP[[I]]
        ## calculate all marker related
        cormat <- cor_col(t(exp), method = "pearson", subset = overlap, ncores = ncores)
        df <- reshape2::melt(cormat)
        df <- df[!is.na(df$value), ]
        df$Genepair <- apply(df[, c(1, 2)], 1, function(x) paste(sort(x), collapse = "~"))
        df <- df[c("Genepair", "value")]
        names(df)[2] <- ds_train[I]
        cordfls[[I]] <- df
        names(cordfls)[I] <- ds_train[I]
        cat("------- Finished\n")
    }
    cordf <- merge.list(cordfls, by = "Genepair", all = T)
    cordf$Type <- apply(cordf[, -1], 1, function(x) {
        x <- x[!is.na(x)]
        numhigh <- sum(abs(x) >= corCutoff)
        consistent <- all(x > 0) | all(x < 0)
        if (numhigh == 0) {
            return("Not-correlated")
        }
        if (numhigh == 1 & consistent) {
            return("Single")
        }
        if (numhigh >= 2 & consistent) {
            return("Consistent")
        }
        return("In-consistent")
    })
    genepair <- as.character(cordf$Genepair[cordf$Type == "Consistent"])
    res <- list(genepair, overlap, EXP, input_markers)
    names(res) <- c("genepair", "markers", "trainScaled","input_markers")
    return(res)
}


pruneNetwork <- function(markers, network = NETWORK){
    nw <- network
    isdef <- identical(nw, NETWORK)
    if(isdef){
        nw$trainScaled <- lapply(CPTAC.6datasets, function(exp) exp <- global.scale(exp)$y)
    }

    if(!is.null(nw$input_markers)){
        if (length(setdiff(markers, nw$input_markers))>0) {
            warnings("Argument 'markers' is not a subset of network's markers, and the overlap will be used!
                     If you don't want to lose information, you'd better to build a new network using buildNetwork() with the latest markers")
        }
        nw$input_markers=intersect(nw$input_markers, markers)
    } else{
        nw$input_markers <- markers
    }


    markers <- intersect(markers, nw$markers)
    gps <- nw$genepair
    gpls <- strsplit(gps, "~")
    gps <- gps[sapply(gpls, function(x) {any(x %in% markers)})]
    nw$markers <- markers
    nw$genepair <- gps
    return(nw)
}

