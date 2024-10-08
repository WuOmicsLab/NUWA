% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NUWA.cibersort.R
\name{NUWA.cibersort}
\alias{NUWA.cibersort}
\title{A portal function running CIBERSORT after NUWAms analysis.}
\usage{
NUWA.cibersort(expr, signature_matrix, cibersortPath, ...)
}
\arguments{
\item{expr}{a numeric matrix of expression profiles for bulk tissue samples, with HUGO gene symbols as rownames and sample identifiers as colnames. Data must be non-logarithm scale.}

\item{signature_matrix}{a signature matrix, which is a numeric expression matrix of markers in cell types of interest, with HUGO gene symbols as rownames and cell type identifiers as colnames. Such as NUWAp26 (a proteomic signature matrix we developed for 26 immune cell types), LM22, LM6, BCIC, TIC or user provided signature matrix.}

\item{cibersortPath}{a string specifying the path of CIBERSORT R script, CIBERSORT is only freely available for academic users, please register on https://cibersort.stanford.edu, and download the CIBERSORT source script.}

\item{...}{additional arguments passed to the NUWAms() function}
}
\value{
The results of each built-in NUWA analysis function, is a list containing an expression matrix with missing markers inferred, two matrices used for recall analysis, and a matrix including cell fractions estimated by the algorithm used.
}
\description{
Run NUWAms (missing markers inference) and CIBERSORT algorithm (deconvolution) with a signature matrix of interest.
}
\examples{
 expr <- CPTAC.6datasets$brca[, 1:5]
 # cibersortPath = "<PATHTO>/CIBERSORT.R"
 res_nuwa <- NUWA.cibersort(expr,
              cibersortPath = cibersortPath,
              signature_matrix = LM22)
 res_nuwa <- NUWA.cibersort(expr,
              cibersortPath = cibersortPath,
              signature_matrix = NUWAp26)
 res_nuwa <- NUWA.cibersort(expr,
              cibersortPath = cibersortPath,
              signature_matrix = my_signature_matrix)

}
