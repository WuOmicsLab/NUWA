# NUWA

## Description

Here, we provide an R package implementing NUWA pipeline for abundance inference of missing function proteins (e.g., cell markers, drug targets) from mass spectrometry-based proteomic profiles, and could lead to improved performance of downstream analyses using proteomic profiles, including deconvolution of immune cell composition, differential expression analysis, etc. 

The major functions are listed below:

(1) `NUWAms`: to infer abundances of missing  proteins based on co-expression networks of individual protein of interest, by leveraging information borrowed from the cohort profiles (training datasets). The default underlying cohort profiles are CPTAC proteomic datasets of six cancer types, which could be replaced by users (e.g., using multiple datasets for a specific cancer type).

(2) `NUWAeDeconv`: a benchmarked ensemble method of three deconvolution algorithm-signature combinations to estimate the relative fractions of six immune cell types. 

(3) A number of portal functions for deconvolution analysis with individual published deconvolution algorithm, following `NUWAms` analysis of proteomic profiles.


## Installation

System requirements: <b>R >= 3.6.1</b>

```R
if(!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github('WuOmicsLab/NUWA')
```
<b>Note 1</b>: if the installation fails due to the `glmnet` package, try manual installation as below before installing `NUWA` package.

Linux users:
`remotes::install_version("glmnet", version = "4.1-1", repos = "https://cran.us.r-project.org")` 

Mac/Windows users:
`install.packages('glmnet', type='binary')`

<b>Note 2</b>: if the installation fails due to the `preprocessCore` package, try manual installation as below before installing `NUWA` package.

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
```

## Usages

The main functions in NUWA package are `NUWAms` and `NUWAeDeconv`. See below for a quick start, while details of each parameters are available in the [manual documentation](https://github.com/WuOmicsLab/NUWA/blob/main/NUWA_1.1.pdf).

### 1) NUWAms

`NUWAms` infers missing values for proteins of interest (termed as marker), with raw proteomic abundance matrix and co-expression networks taken as input. 

(**a**) Run `NUWAms` using the default co-expression networks for 10487 proteins (`NETWORK_LIST.10487markers`), which were constructed using CPTAC proteomic datasets of six cancer types (`CPTAC.6datasets`) including breast cancer, clear cell renal cell carcinoma, colon cancer, endometrial cancer, gastric cancer and lung cancer. 

```R
library(NUWA)  # load NUWA package
res_nuwams <- NUWAms(expr = raw_expr, # Provided raw proteomic abundance matrix for inference of missing values
                    markers = my.markers, # Provided proteins of interest
                    )
```
`res_nuwams` is a list including an expression matrix after abundance inference of missing markers, and additional matrices used to evaluate the inference accuracy. 

<b>Note</b>: NUWAmsonly infers abundance for markers having a co-expression network. 

(**b**) Run `NUWAms` uing the marker co-expression networks, bulit using user provided training datasets (e.g. multiple datasets for a specific cancer type). An additional step is then needed to  build co-expression networks, by running function `buildNetwork`.

```R
my.network <- buildNetwork(trainsets = list_trainsets, # Provided list containing user provided multiple trainsets 
                  markers = my.markers)

res_nuwams <- NUWAms(expr = raw_expr,
                network = my.network, # The customized markers' co-expression networks built in the previous step
                markers = my.markers)
```

### 2) NUWAeDeconv

`NUWAeDeconv` takes the `NUWAms` inferred protein abundance matrix, and the path of CIBERSORT R script as input. 

<b>Note</b>: CIBERSORT is freely available for academic users, but users are needed to register on [its website](https://cibersort.stanford.edu) to download the source script.

```R
cibersortPath = "<PATHTO>/CIBERSORT.R" # The local file path to CIBERSORT R source code

res_deconv <- NUWAeDeconv(expr = res_nuwams$finalExpr, cibersortPath = cibersortPath)
```

`res_deconv` includes matrices for immune cell fractions estimated by `NUWAeDeconv`,  original predictions and updated ones (with cell types merged)  by CIBERSORT-LM22, CIBERSORT-LM6 and EPIC-BCIC, respectively. 


### 3) Other deconvolution approaches

Additional convenient portal functions are provided for deconvolution analysis using individual published deconvolution algorithm (including CIBERSORT, EPIC, MCPcounter and xCell), following `NUWAms` analysis of proteomic profiles.

```R
# run NUWAms and EPIC
res_nuwa <- NUWA.EPIC(expr = raw_expr, signature_matrix = BCIC)

# run NUWAms and xCell algorithm
res_nuwa <- NUWA.xcell(expr = raw_expr)

# run NUWAms and MCPcounter algorithm
res_nuwa <- NUWA.mcpcounter(expr = raw_expr)
```
`res_nuwa` is a list, including an expression matrix after abundance inference of missing markers using NUWAms modelling, and a matrix with immune cell fractions estimated by the selected algorithm.

## License
NUWA is free for academic users of non-commercial research. Commercial use of NUWA requires a license (contact Dr. Jianmin Wu by wujm@bjmu.edu.cn for details). 

## Metaphor of the package name 'NUWA'
In Chinese mythology, [Nuwa](https://mythopedia.com/topics/nuwa) is considered to be the first being and has a famous story for saving humanity by mending a hole in the sky.

## Contact information
Lihua Cao (lihuacao@bjcancer.org), Yuhao Xie (xieyuhao@pku.edu.cn), and Jianmin Wu (wujm@bjmu.edu.cn).
