# NUWA

## Description

An R package implementing NUWA pipeline for abundance inference of missing cell markers in mass spectrometry-based proteomic profiles, to enable accurate deconvolution of relative immune cell fractions. It consists of three major modules: (1) `NUWAms`, to infer abundances of missing cell markers based on co-expression networks of individual markers. For users' convenience, the package includes built-in cancer co-expression networks for markers of <b>NUWAp26</b> (a proteomic signature matrix we developed for 26 immune cell types), and a set of previously published markers (signature genes), including BCIC, LM22, LM6, MCPcounter and xCell. (2) `NUWAeDeconv`, a benchmarked ensemble method of three deconvolution algorithm-signature combinations to estimate the relative fractions of six immune cell types. (3) A number of portal functions for deconvolution analysis with individual published deconvolution algorithm, following `NUWAms` analysis of proteomic profiles.


## Installation

System requirements: R >= 3.6.1

```R
if(!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github('WuOmicsLab/NUWA')
```
Note: if the installation fails due to `glmnet` package, try manual installation by `remotes::install_version("glmnet", version = "4.1-1", repos = "https://cran.us.r-project.org")` before installing `NUWA` package. 

## Usages

The main functions in NUWA package are `NUWAms` and `NUWAeDeconv`. See below for a quick start, while details of each parameters are available in the [manual documentation](https://github.com/WuOmicsLab/NUWA/blob/main/man/NUWA_0.1.0.pdf).

### 1) NUWAms

`NUWAms` takes a protein abundance matrix and a co-expression networks (optional) as input. Each column of protein abundance matrix represents a sample and each row represents a protein.  

(**a**) Run `NUWAms` using the default co-expression network for LM22, LM6 and BCIC signatures:

```R
library(NUWA)  ## load NUWA package
res_nuwams = NUWAms(expr = raw_expr, network = NULL)
```
"res_nuwams" is a list including an expression matrix after abundance inference of missing markers using NUWAms modelling, and additional matrices used to evaluate the inference accuracy.

(**b**) Run `NUWAms` with a customized co-expression network, which is built by function `buildNetwork` using user provided training datasets and marker genes:

```R
my.network <- buildNetwork(trainsets = cptacDatasets, markers = my.markers)
res_nuwams <- NUWAms(expr = raw_expr, network = my.network)
```

### 2) NUWAeDeconv

`NUWAeDeconv` takes the `NUWAms` inferred protein abundance matrix, and the path of CIBERSORT R script as input. Note: CIBERSORT is  freely available for academic users, but request users  to register on its website [https://cibersort.stanford.edu](https://cibersort.stanford.edu) to download the  source script.

```R
# Provide the file path to CIBERSORT R source code
cibersortPath = "<PATHTO>/CIBERSORT.R"
res_deconv <- NUWAeDeconv(expr = res_nuwams$finalExpr, cibersortPath = cibersortPath)
```

"res_deconv" includes matrices for immune cell fractions estimated by `NUWAeDeconv`,  original predictions and updated ones (with cell types merged)  by CIBERSORT-LM22, CIBERSORT-LM6 and EPIC-BCIC, respectively. 


### 3) Other deconvolution approaches

Additionally, convenient portal functions are provided for deconvolution analysis using individual published deconvolution algorithm (including CIBERSORT, EPIC, MCPcounter and xCell), following `NUWAms` analysis of proteomic profiles.

```R
# run NUWAms and EPIC
res_nuwa <- NUWA.EPIC(expr = raw_expr, signature_matrix = BCIC)

# run NUWAms and xCell algorithm
res_nuwa <- NUWA.xcell(expr = raw_expr)

# run NUWAms and MCPcounter algorithm
res_nuwa <- NUWA.mcpcounter(expr = raw_expr)
```
"res_nuwa" is a list, including an expression matrix after abundance inference of missing markers using NUWAms modelling, and a matrix with immune cell fractions estimated by the selected algorithm.

## License

NUWA is free for academic users of non-commercial purposes. Commercial use of NUWA requires a license. If NUWA package was used for your analysis, please cite our package and the used deconvolution algorithm(s).


## References
| Name | license | PMID | Citation |
| :- | :- | :- |:- |
| MCPcounter | free | 27765066 | Becht, E. et al. Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biol 17, 218 (2016). |
| xCell | free | 29141660 | Aran, D., Hu, Z. & Butte, A.J. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol 18, 220 (2017). |
| CIBERSORT | free for non-commercial use | 25822800 | Newman, A.M. et al. Robust enumeration of cell subsets from tissue expression profiles. Nat Methods 12, 453-457 (2015). |
| EPIC | free | 29130882 | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D.E. & Gfeller, D. Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. Elife 6 (2017). |

## Contact information
Lihua Cao (lihuacao@bjcancer.org), Yuhao Xie (xieyuhao@pku.edu.cn), Jianmin Wu (wujm@bjmu.edu.cn).
