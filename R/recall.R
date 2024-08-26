#' Evaluate the inference accuracy of NUWAms by recall analysis.
#'
#' This function computes recall value for each marker and sample based on the inference results of NUWAms function or NUWA built-in analysis functions. 
#' Recall of a marker was determined as the fraction of a null distribution of similarity less than the observed similarity (Spearman rank correlation)
#' between the inferred and measured abundances of this marker in all samples of a given dataset.
#' Recalls cannot be computed if the number of markers with both NUWA-ms inference and quantification in the inferred dataset is less than 10.
#' Recall of a sample was performed in a similar way. 
#'
#' @param x a list, the output of NUWAms().
#' @param corMethod a character specifying the type of correlation coefficient to be used: "pearson" or "spearman" (default).
#'
#' @return A list containing: \describe{
#'  \item{\code{recallTable}}{a data frame recording the recall values of each marker and each sample.}
#'  \item{\code{simNull}}{a list of length 2, consisting two correlation matrices which were used to generate null distributions of similarity (SIMnull) at marker and sample level, respectively. See the Methods section of the NUWA manuscript for more detail.}
#'  \item{\code{corMethod}}{a character, the type of correlation method used.}
#'  }
#' @export
#'
#' @examples
#' expr <- CPTAC.6datasets$brca[, 1:20]
#' myInfer <- NUWAms(expr)
#' recallRes <- recall(myInfer)

recall <- function(x, corMethod = c("spearman", "pearson")[1]){
    truth.pred=x$predVsTruth
    pred=truth.pred[['pred']]
    if(min(dim(pred))<=3) {
        stop("Recall can not be computed when the number of samples or genes is less than or equal to 3.")
    }
    truth=truth.pred[['truth']]

    geneIndex=1:nrow(truth)
    names(geneIndex)=rownames(truth)
    sampleIndex=1:ncol(truth)
    names(sampleIndex)=colnames(truth)
    rmHighCor <- T

    if(rmHighCor){
        truth.cor.gene <- sapply(geneIndex, function(x){
            xx <- truth[x,]
            res <- sapply(geneIndex, function(y){
                if (x>=y) return(0) else{
                    yy <- truth[y,]
                    return(cor(xx,yy,use='p',method=corMethod))
                }
            })
        })
        truth.cor.gene <- truth.cor.gene + t(truth.cor.gene)
    }
    simNull.cor.gene=sapply(geneIndex, function(x){
        xx=truth[x,]
        res=sapply(geneIndex,function(y){
            if(x==y) {
                return(NA)
            }else if (rmHighCor && !is.na(truth.cor.gene[x,y]) && truth.cor.gene[x,y]>0.5){
                return(NA)
            } else{
                yy=pred[y,]
                return(cor(xx,yy,use='p',method=corMethod))
            }

        })
    })
    simNull.cor.sample=sapply(sampleIndex, function(x){
        xx=truth[,x]
        res=sapply(sampleIndex,function(y){
            if(x==y) return(NA) else{
                yy=pred[,y]
                return(cor(xx,yy,use='p',method=corMethod))
            }

        })
    })

    simNull.list=list(Gene=simNull.cor.gene,
                      Sample=simNull.cor.sample)

    df.recall.ls=lapply(c('Gene','Sample'), function(level){
        if(level == 'Sample'){
            pred0=t(pred)
            truth0=t(truth)
        }else{
            pred0=pred
            truth0=truth
        }
        simNull=simNull.list[[level]]
        simNull=as.vector(simNull)
        simNull=simNull[!is.na(simNull)]
        gs=rownames(pred0)
        names(gs)=gs
        dfls=lapply(gs,function(g){
            x=pred0[g,]
            y=truth0[g,]
            cor=cor(x,y,use='p',method = corMethod)
            pval=sum(simNull>cor)/length(simNull)
            recall=1-pval
            res=data.frame(g,cor,pval,recall,level)
            names(res)=c('sample','cor','pval','recall','level')
            return(res)
        })
        df=do.call(rbind,dfls)
        return(df)
    })

    df.recall=do.call(rbind, df.recall.ls)
    rownames(df.recall)=NULL


    finalRes=list(
        recallTable=df.recall,
        simNull=simNull.list,
        corMethod=corMethod
    )
    recallsum(df.recall)
    return(finalRes)
}
