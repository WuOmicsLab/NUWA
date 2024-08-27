#' Infer proteomic abundance using the given co-expression networks of individual marker, for user provided marker proteins.
#'
#' This function is to infer abundances of missing proteins based on co-expression networks of individual protein of interest (termed as marker), by leveraging information borrowed from the cohort profiles (training datasets). The default underlying cohort profiles are CPTAC proteomic datasets of six cancer types (breast cancer, clear cell renal cell carcinoma, colon adenocarcinoma, endometrial carcinoma, gastric cancer, lung adenocarcinoma), which could be replaced by users (e.g., using multiple datasets for a specific cancer type). It may take a few hours, if more than 50 samples were included in the analysis.
#'
#' @param expr a numeric matrix of expression profiles for bulk tissue samples, with HUGO gene symbols as rownames and sample identifiers as colnames.
#' 
#' @param network a list, consisting of the co-expression networks built by function buildNetwork(). A list containing built-in networks for 10487 proteins (NETWORK_LIST.10487markers) will be used, when Default (NULL) is applied. It was constructed using CPTAC proteomic datasets of six cancer types. Please note, we only infer abundance for markers existing in the provided network. 
#' 
#' @param markers a character vector of interesting proteins (termed as marker) to infer. If NULL (default) provided, a pre-defined list MARKERS_Immune.1114 (1114 immune marker proteins collected from public studies) will be used. Another two pre-defined lists, including MARKER_FDA_Drug_Targets.812 (FDA-approved drug targets from HPA and Drugbank) and MARKERS_Immune.NUWAp26 (631 immune cell markers from proteomic signature matrix NUWAp26 we developed for 26 immune cell types) were also provided for convenience.
#' 
#' @param direction a character, indicating the mode used for feature searching in stepwise regression analysis, one of "both", "backward" or "forward". Default is "both".
#' @param lasso_step_cutoff a positive integer, specifying the minimal number of variables needed to run LASSO regression analysis. If the number is less than "lasso_step_cutoff", stepwise regression models will be constructed. Default is 10.
#' 
#' @param preprocess logical. If TURE, expression data is preprocessed before markers inference. Default is TRUE. See preprocess() function and the Methods section of the NUWA manuscript for more details.
#' @param ncores a positive integer, indicating the number of cores used by this function. If the operating system is windows, then only one core will be used.
#' @param lambda a character, indicating which value of lambda will be used in the LASSO analysis. One of "lambda.min" or "lambda.1se". "lambda.min" gives lambda with minimal cross-validation errors, and "lambda.1se" gives the largest value of lambda such that the error is within 1 standard error of the minimal. Default is "lambda.1se".
#' 
#' @param quantification_method The quantification method used in MS-proteomic profiling, "TMT/iTRAQ ratio" or "Label-free intensity". Default is "TMT/iTRAQ ratio".
#' @param logbase The log base of the expression matrix if the quantification values have been log-transformed, one of "No transformation" or "Log2". Default is "No transformation".

#' @return A list containing:\describe{
#'  \item{\code{finalExpr}}{a numeric matrix, the final full dataset of expression with missing markers are inferred.}
#'  \item{\code{predVsTruth}}{a list with elements comprising the prediction and truth expression matrices of quantified markers, which will be used for the following recall analysis.}
#'  \item{\code{inferenceMat}}{a numeric matrix, a subset of the full dataset "finalExpr", a expression matrix only including markers inferred by NUWAms() function.}
#'  \item{\code{runtime}}{a data frame, consisting running time (minutes) used for model building, markers inferring and total time.}
#' }
#'
#' @export
#'
#' @examples
#' expr <- CPTAC.6datasets$brca[, 1:5]
#' res <- NUWAms(expr)
#' res <- NUWAms(expr, markers = MARKER_FDA_Drug_Targets.812)
#' 

NUWAms <- function(expr, network = NULL, markers = NULL,
                   direction=c("both", "backward", "forward")[1],
                   lasso_step_cutoff=10, preprocess = T, ncores = 16,
                   lambda = c("lambda.1se", "lambda.min")[1],
                   quantification_method=c("TMT/iTRAQ ratio", "Label-free intensity")[1],
                   logbase = c("No transformation", "Log2")[1]
                   ){

    ##
    simplifyModel <- function(mod) {
        if(is.null(mod)) return(NULL)
        cla <- class(mod)[1]
        if(cla == "lm") {
            coef <- coef(mod)
            res <- list(class = cla, coef = coef)
        } else if (cla == "cv.glmnet") {
            coef <- coef(mod,s = c(mod$lambda.1se, mod$lambda.min))
            res <- list(class = cla, coef = coef)
        } else {
            warning("Unrecognized model!")
            res <- NULL
        }
    }

    ##
    mypredict <- function(mod, x, lambda = c("lambda.1se", "lambda.min")[1]) {
        cla <- mod$class
        coef0 <- mod$coef
        if(cla ==  "lm") {
            x0 <- data.matrix(cbind(1, x[names(coef0)[-1]]))
        } else if (cla == "cv.glmnet") {
            library("Matrix")
            x0 <- data.matrix(cbind(1, x[rownames(coef0)[-1]]))
            coef0 <- if (lambda == "lambda.1se") coef0[, 1, drop = F] else coef0[, 2, drop = F]
        }
        res <- x0 %*% coef0
        res <- res[, 1]
        return(res)
    }

    # mar=Marker[,1]
    # set.seed(1001) 
    
    # non-log trans
    if(any(expr < 0, na.rm = T) | logbase == "Log2") {
        expr <- 2^expr
    }

    # if quantification_method is "label-free-DIA", must set preprocess = TRUE
    if (quantification_method == "Label-free intensity") { preprocess = TRUE } 

    prep=preprocess
    method <- c('lasso_step','step','RF')[1]
    lazymode <- F
    t0 <- Sys.time()
    is.window <- (get_os() == "windows")
    if (is.window) ncores <- 1 else ncores <- min(ncores, parallel::detectCores()-1)
    expr=data.matrix(expr)
    if(prep) {
        cat("Data preprocessing \n")
        expr=preprocess(expr, thre = 0, quantification_method = quantification_method)
        cat("Finished!\n\n")
    }
    expr_gsls <- global.scale(expr)
    expr <- expr_gsls$y
    reverse_fun <- expr_gsls$reverseFun
    # c("corr", "markers", "trainScaled")
    if(is.null(markers)) markers <- MARKERS_Immune.1114
    if (is.null(network)) network <- NETWORK_LIST.10487markers
    network <- pruneNetwork(markers, network)


    markers <- unique(network$markers)
    EXP <- network$trainScaled
    # cordf <- network$corr
    # gps <- cordf$Genepair[cordf$Type == 'Consistent']
    gps <- network$genepair
    gpls <- strsplit(gps, "~")
    gpls <- gpls[sapply(gpls, length) == 2]
    gpmat <- do.call(rbind, gpls)
    gpmat <- rbind(gpmat, gpmat[, c(2, 1)])
    cor_genes <- unique(gpmat[, 1])
    cor_marker <- intersect(cor_genes, markers)
    df.genepair <- data.frame(gpmat)
    names(df.genepair) <- c("x", "y")
    df.genepair <- df.genepair[df.genepair$y %in% cor_marker, ]
    sams=colnames(expr)
    names(sams)=sams
    rescale <- T
    expr1=expr
    t1 <- Sys.time(); print(t1 - t0)

    # part 1: constructing models.------------------------------
    cat("Building regression models\n")
    if (is.window) cl <- parallel::makePSOCKcluster(ncores) else cl <- parallel::makeForkCluster(ncores)
    # parallel::clusterSetRNGStream(cl, 2019)
    file <- paste0(getwd(),"/model_", ncol(expr), "_", ncores,".rda")
    # lazymode <- F
    if(file.exists(file) && lazymode) {load(file)} else{
        # tryCatch({}, warning = function(w) print(w))

        modls=lapply(EXP, function(df){
            df1=data.frame(t(df))
            gene.tr=colnames(df1)
            gene.tr.y=intersect(gene.tr,cor_marker)
            names(gene.tr.y)=gene.tr.y
            if(method %in% c('step','RF')){
                gene.tr.x.ls=lapply(sams,function(s){
                    gx=unique(rownames(expr1)[!is.na(expr1[,s])])
                    gx=intersect(gx,gene.tr)
                })
                mls=parallel::parLapply(cl, gene.tr.y,function(gy){
                    gx=df.genepair$x[df.genepair$y==gy]
                    formula.x=sapply(gene.tr.x.ls,function(gene.tr.x){
                        gene.x=intersect(gene.tr.x,gx)
                        gene.x=sort(unique(gene.x))
                        matx=data.matrix(df1[!is.na(df1[,gy]),gene.x])
                        gene.x.index=apply(matx,2,function(x){all(!is.na(x))})
                        gene.x=gene.x[gene.x.index]
                        if(length(gene.x)==0) fx=NA else fx=paste0(gene.x,collapse = '+')
                        return(fx)
                    })
                    formula.x.uniq=unique(formula.x[!is.na(formula.x)])
                    model.index=match(formula.x, formula.x.uniq)
                    model.ls=lapply(formula.x.uniq,function(fx){
                        formula=as.formula(paste0(gy,'~',fx))
                        if(method=='step'){
                            res=tryCatch({
                                lm.mod=lm(formula,df1)
                                step.mod=step(lm.mod,direction = direction, trace = F)
                            },error=function(e){return(NULL)})
                        }else if(method=='RF'){
                            res=randomForest::randomForest(formula,df1,na.action=na.omit)
                        }
                        return(simplifyModel(res))
                    })
                    res1=list(modelIndex=model.index,modelList=model.ls)
                    return(res1)
                })
                # tryCatch({}, warning = function(w) print(w))
                # tryCatch({}, warning = function(w) print(w))
            }else if(method=='lasso_step'){
                gene.tr.x.ls=lapply(sams,function(s){
                    gx=unique(rownames(expr1)[!is.na(expr1[,s])])
                    gx=intersect(gx,gene.tr)
                })
                # tryCatch({}, warning = function(w) print(w))
                # tryCatch({}, warning = function(w) print(w))
                mls=parallel::parLapply(cl, gene.tr.y, function(gy){
                    gx=df.genepair$x[df.genepair$y==gy]
                    formula.x=sapply(gene.tr.x.ls,function(gene.tr.x){
                        gene.x=intersect(gene.tr.x,gx)
                        gene.x=sort(unique(gene.x))
                        matx=data.matrix(df1[!is.na(df1[,gy]),gene.x])
                        gene.x.index=apply(matx,2,function(x){all(!is.na(x))})
                        gene.x=gene.x[gene.x.index]
                        if(length(gene.x)==0) fx=NA else fx=paste0(gene.x,collapse = '+')
                        return(fx)
                    })
                    formula.x.uniq=unique(formula.x[!is.na(formula.x)])
                    model.index=match(formula.x, formula.x.uniq)
                    model.ls=lapply(formula.x.uniq,function(fx){
                        gx=strsplit(fx,'\\+')[[1]]
                        if(length(gx)<=lasso_step_cutoff){
                            # print('step lm method')
                            formula=as.formula(paste0(gy,'~',fx))
                            res=tryCatch({
                                lm.mod=lm(formula,df1)
                                step.mod=step(lm.mod,direction = direction, trace = F)
                            },error=function(e){print(e);return(NULL)})
                        }else{
                            # print('lasso method')
                            set.seed(2019)
                            matx=data.matrix(df1[,gx])
                            maty=df1[,gy]
                            mat=cbind(matx,maty)
                            na.index=apply(mat,1,function(x) any(is.na(x)))
                            matx=matx[!na.index,]
                            maty=maty[!na.index]
                            if(length(maty)>=30) nfolds=10 else nfolds=floor(length(maty)/3)
                            if(nfolds<=3) res=NULL else{
                                res <- glmnet::cv.glmnet(x=matx, y=maty, alpha = 1, family='gaussian', nfolds=nfolds)
                            }
                        }
                        return(simplifyModel(res))

                    })
                    res1=list(modelIndex=model.index,modelList=model.ls)
                    return(res1)
                })
            }
            return(mls)
        })
        # parallel::stopCluster(cl)
        # nn8=modls11[[1]][[4]][[2]][[2]][[2]][,1]
        # nn8[nn8!=0]
        # save(list=c('modls'),file = file)
    }
    t2 <- Sys.time(); print(t2 - t1)
    cat("Finished!\n\n")
    model_time <- difftime(t2, t0, units = "mins")

    ## part 2: filling ---------------------------------
    cat("Inferring Markers \n")
    modls=lapply(modls,function(mls) {
        mls[sapply(mls,function(x) identical(x,NULL))]=NULL
        return(mls)})
    gene.y=unique(unlist(lapply(modls,names)))
    names(gene.y)=gene.y
    # print(str(gene.y))
    expr2=data.frame(t(expr1))

    if(method == "RF") {
        library(randomForest)
    }
    expr.fill.ls=lapply(modls, function(mls){
        # expr.fill.ls=lapply(modls, function(mls){
        y=parallel::parSapply(cl, gene.y,function(g){
            yy=rep(NA,nrow(expr2))
            if(g %in% names(mls)){
                gene.model.all=mls[[g]]
                model.index=gene.model.all$modelIndex
                gene.model.ls=gene.model.all$modelList
                if(length(gene.model.ls)==0) return(yy)
                for(i in 1:length(gene.model.ls)){
                    gene.model=gene.model.ls[[i]]
                    if(!identical(gene.model,NULL)){
                        expr3=expr2[model.index==i & !is.na(model.index),]
                        if(class(gene.model)[1] == "randomForest.formula"){
                            yy[model.index==i & !is.na(model.index)]=predict(gene.model,expr3)
                        }else {
                            yy[model.index==i & !is.na(model.index)]=mypredict(gene.model,expr3,lambda)
                        }
                    }
                }
            }
            return(yy)
        })
        if(is.vector(y)) y <- t(y)
        rownames(y)=rownames(expr2)
        return(t(y))
    })
    parallel::stopCluster(cl)
    #   parallel::clusterExport(cl,varlist = setdiff(ls(), "cl"), envir = environment())
    # tryCatch({}, warning = function(w) {cat("Filling error! outside \n");print(w)})
    t3 <- Sys.time(); print(t3 - t2)
    cat("Finished!\n\n")
    # save(file=filled.file, list = c("expr.fill.ls"))

    # Part 3. remove outliers for each sample -------------------
    cat("Removing outliers\n")
    expr.fill.ls <- lapply(
        expr.fill.ls,
        function(matr) {
            for (i in 1:ncol(expr)) {
                x <- expr[,i]
                x <- x[!is.na(x)]
                qnt <- quantile(x, probs= c(.25,.75))
                H <- 1.5*IQR(x)
                min0 <- 0
                max0 <- max((qnt[2]+H), quantile(x,0.99))
                y <- matr[,i]
                y[y>max0|y<min0] <- NA
                matr[,i] <- y
            }
            return(matr)
        }
    )
    cat("Finished!\n\n")
    # part 3.2: integrate multiple outputs into one ----------------
    expr.fill=abind::abind(expr.fill.ls,along = 3)
    expr.fill=apply(expr.fill, c(1,2), mean,na.rm=T)
    expr.fill=expr.fill[!apply(is.na(expr.fill),1,all),, drop = F]


    # part4. rescale by sample ----------------------
    genes0 <- intersect(rownames(expr.fill), rownames(expr))
    if (length(genes0)<20) {
        warning("The overlaped markers between truth and predictions are ",
            "fewer than 20, so rescaling process is skipped to maintain ",
            "numeric stablity")
        rescale <- FALSE
    }
    if(rescale) {
        cat("Rescaling\n")
        mat.pred <- expr.fill[genes0, ,drop=F]
        mat.truth <- expr[genes0, ,drop=F]
        for(s in colnames(expr)) {
            # print(s)
            pred <- mat.pred[,s]
            truth <- mat.truth[,s]
            mean.p <- mean(pred, na.rm= T)
            mean.t <- mean(truth, na.rm=T)
            std.p <- sd(pred, na.rm = T)
            std.t <- sd(truth, na.rm = T)
            expr.fill[,s] <- (expr.fill[,s]-mean.p)*(std.t/std.p)+mean.t
            rescale <- expr.fill[genes0, s]
        }
        cat("Finished!\n\n")

    }

    # part 5. expr.final integration-----------------------
    cat("Integrating\n")
    getValueFromMat=function(matrix,rowname,colname){
        if(rowname %in% rownames(matrix) & colname %in% colnames(matrix)){
            return(matrix[rowname,colname])
        }else{
            return(NA)
        }
    }
    genes.all=sort(unique(union(rownames(expr),rownames(expr.fill))))
    names(genes.all)=genes.all
    expr.final=sapply(sams,function(s){
        y=sapply(genes.all,function(g){
            value.truth=getValueFromMat(expr,g,s)
            value.pred=getValueFromMat(expr.fill,g,s)
            if(!is.na(value.truth)){
                yy=value.truth
            }else{
                if(is.na(value.pred)) yy=NA else yy=value.pred
            }
            return(yy)
        })
    })
    expr.final=expr.final[!apply(is.na(expr.final),1,all),, drop = F]
    cat("Finished!\n\n")


    # part 6. pred vs truth-----------------------------
    gene.common=intersect(rownames(expr.fill),rownames(expr))
    truth.pred=list(truth=expr[gene.common,, drop = F],pred=expr.fill[gene.common,, drop = F])
    t4 <- Sys.time(); print(t4 - t3)
    cat("total:"); print(t4-t0)
    infer_time <- difftime(t4, t2, units = "mins")
    total_time <- difftime(t4, t0, units = "mins")
    nsams <- ncol(expr)
    runtime_mins <- data.frame(
        nsams = nsams,
        ncore = ncores,
        model_time = as.numeric(model_time),
        infer_time = as.numeric(infer_time),
        total_time = as.numeric(total_time)
    )
    finalExpr = reverse_fun(expr.final)
    finalExpr[finalExpr < 0] <- 0
    finalRes <- list(
        finalExpr = finalExpr,
        predVsTruth = lapply(truth.pred, reverse_fun),
        inferenceMat = reverse_fun(expr.fill),
        runtime = runtime_mins
    )
    return(finalRes)
}

