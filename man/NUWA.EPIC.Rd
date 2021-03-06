% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NUWA.EPIC.R
\name{NUWA.EPIC}
\alias{NUWA.EPIC}
\title{A portal function running EPIC after NUWAms analysis.}
\usage{
NUWA.EPIC(expr, signature_matrix, ...)
}
\arguments{
\item{expr}{see the same argument in NUWA.cibersort.}

\item{signature_matrix}{see the same argument in NUWA.cibersort.}

\item{...}{additional arguments passed to the NUWAms() function}
}
\value{
see NUWA.cibersort.
}
\description{
Run NUWAms and EPIC algorithm with a signature matrix of interest.
}
\examples{
expr <- cptacDatasets$brca[, 1:5]
res_nuwa <- NUWA.EPIC(expr, signature_matrix = BCIC)
res_nuwa <- NUWA.EPIC(expr, signature_matrix = TIC)
res_nuwa <- NUWA.EPIC(expr, signature_matrix = my_signature_matrix)
}
