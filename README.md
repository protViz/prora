# fgcz.gsea.ora

An R package for running enrichment analysis.  

A plethora of R packages exist on CRAN and Bioconductor to perform over-representation 
analysis (ORA) and gene set enrichment analysis (GSEA). However, consistency in the 
underlying nomenclature for specific analyses and user friendly implementation is 
still lacking. `fgcz.gsea.ora` aims at unifying ID mapping and enrichment analysis 
in a syntactically coherent and intuitive way, while ensuring reproducibility of 
results. 

`fgcz.gsea.ora` primarily consists of wrapper functions around the 
`r CRANpkg("sigora")` and `r CRANpkg("WebGestaltR")` packages from CRAN and 
`r CRANpkg("rmarkdown")` based reports for visualisation and contextualisation 
of analysis results. 

# 1. Future plans

1) Integration of the R packages `r Biocpkg(topGo)` and  `r CRANpkg("enrichr")`
2) Providing executable R script files, which can be run on windows or linux.

# 2. System requirements

## 2.1 Currently the software have been tested on the following systems


| platform | platform version  | R version  |
--- | --- | --- |
| Windows  | 10 x64          | R-3.6.2      |


# 3. Installation guide:

# 3.1 Dependancies

Run the following lines in R to make sure all required ackages are installed

```
list_of_packages <- c("AnnotationDbi", "BiocStyle", "dplyr", "DT", "ggplot2", "GO.db", "httr", "magritrr", "org.Hs.eg.db", "reactome.db", "readr", "rlang", "S4Vectors", "sigora", "slam", "tibble", "tidyr", "tidyverse", "UpSetR")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

## 3.2 **Windows:**

Run the following in R:

```
install.packages("remotes")
library(remotes)
remotes::install_github("protViz/fgcz.gsea.ora", build_vignettes = FALSE)
```

For installing package together with vignettes run:
```
remotes::install_github("protViz/fgcz.gsea.ora", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```


# 4. Vignettes

See fgcz.gsea.ora vignette:

```
vignette("fgcz.gsea.ora")
```

Vignette can also be found here:

https://rdrr.io/github/protViz/fgcz.gsea.ora/f/vignettes/vignette.Rmd


Tools to look at.

[https://agotool.org/]
[https://www.biorxiv.org/content/10.1101/731596v2]

