# drugMap
### Pipeline description

This repository contains scripts to run the analysis pipeline described in the manuscript:

"Integrated network pharmacology approach for drug combination discovery: a multi-cancer case study" by Antonio Federico, Michele Fratello, Giovanni Scala, Lena Möbus, Alisa Pavel, Giusy del Giudice, Michele Ceccarelli, Valerio Costa, Alfredo Ciccodicola, Vittorio Fortino, Angela Serra and Dario Greco (submitted).

The pipeline allows to perform advanced network modelling operations in order to reproduce the results reported in the aforementioned manuscript.

### Dependency libraries 

nettools R package available at: https://github.com/filosi/nettools
Installation from GitHub:

```R
install.packages("devtools")
library(devtools)
install_github("filosi/nettools")
```


#### Download L1000 data from GEO
In order to run the analyses L1000 data and metadata are needed. The files are available on GEO and can be downloaded from their ftp server.
The following command can be used in R to download the data.

```R
remotes::install_github("skgrange/threadr")
library(threadr)

download_ftp_file(file_remote = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz",
                  file_local = "/home/MOKA/Angela/angela/drugMap/data/GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx.gz")

download_ftp_file(file_remote = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                  file_local = "/home/MOKA/Angela/angela/drugMap/data/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")

download_ftp_file(file_remote = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz",
                  file_local = "/home/MOKA/Angela/angela/drugMap/data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz")

download_ftp_file(file_remote = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz",
                  file_local = "/home/MOKA/Angela/angela/drugMap/data/GSE92742_Broad_LINCS_sig_info.txt.gz")

```

