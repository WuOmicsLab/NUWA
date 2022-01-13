# NUWA

## Description

An R package implementing NUWA pipeline to infer abundances of missing cell markers and decipher relative immune cell fractions using mass spectrometry-based proteomic profiles. It consists of two major functions: (1) `NUWAms`,  to infer  abundances of missing cell markers using the given co-expression networks of individual marker. For users' convenience, the package includes built-in cancer co-expression networks for markers of BCIC, LM22, LM6, MCPcounter and xCell signatures/marker lists. (2) `NUWAeDeconv`, an ensemble method of three deconvolution algorithm-signature combinations to estimate the relative fractions of six immune cell types. 

## Installation

System requirements: R >= 3.6 

```R
install.packages("devtools")
devtools::install_github('Wulab-CCB/NUWA')
```

## Usages

The main functions in NUWA package are `NUWAms`  and `NUWAeDeconv`. See below for a quick start, while details of each parameters is available at the document file.

### 1) NUWAms

NUWAms takes a protein abundance matrix and a co-expression networks (optional) as input. Each column of protein abundance matrix represents a sample and each row represents a protein.  

(**a**) Run `NUWAms` using the default co-expression network for LM22, LM6 and BCIC signatures:

```R
library(NUWA)  ## load NUWA package
res_nuwams = NUWAms(expr = raw_expr, network = NULL)
```
"res_nuwams" is a list including an expression matrix after abundance inference of missing markers using NUWAms modelling, and additional matrices used to evaluate the inference accuracy.

(**b**) Run `NUWAms` with a customized co-expression network, which is built by function `buildNetwork` using user provided training datasets and marker genes:

```R
my.network <- buildNetwork(trainsets = cptacDatasets, markers = my.markers)
res_nuwams = NUWAms(expr = nuwa_raw_expr, network = my.network)
```

### 2) NUWAeDeconv

`NUWAeDeconv` takes the `NUWAms` inferred protein abundance matrix, and the path of CIBERSORT R script as input. Note: CIBERSORT is  freely available for academic users, but request users  to register on its website [https://cibersort.stanford.edu](https://cibersort.stanford.edu) to download the  source script.

```R
# Provide the file path to CIBERSORT R source code
cibersortPath = "<PATHTO>/CIBERSORT.R"
res_deconv = NUWAdeconv(m.exp = res_nuwams$finalExpr, cibersortPath = cibersortPath)
```

"res_deconv" includes matrices for immune cell fractions estimated by `NUWAeDeconv`,  originial predictions and updated ones (with cell types merged)  by CIBERSORT-LM22, CIBERSORT-LM6 and EPIC-BCIC, respectively. 


### 3) other deconvolution approaches

Beside the ensemble deconvolution method `NUWAeDeconv`, the package also supports analysis using individual algorithm (including EPIC, Cibersort, MCPcounter, xCell) using the built-in co-expression networks for markers of these interested algorithms.

```R
# run NUWAms and EPIC
res_nuwa <- NUWA.EPIC.BCIC(expr = raw_expr, protein = TRUE)

# run NUWAms and xCell algorithm
res_nuwa <- NUWA.xcell.64(expr = raw_expr)

# run NUWAms and MCPcounter algorithm
res_nuwa <- NUWA.mcpcounter(expr = raw_expr)
```
"res_nuwa" is a list, including an expression matrix after abundance inference of missing markers using NUWAms modelling, and a matrix with immune cell fractions estimated by the selected algorithm.

## License

NUWA is free for academic users of non-commercial purposes. Commercial use of NUWA requires a license.

## References
