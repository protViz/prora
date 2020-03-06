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

# Future plans

1) Integration of the R packages `r Biocpkg(topGo)` and  `r CRANpkg("enrichr")`
2) Providing executable R script files, which can be run on windows or linux.

# Installation guide:

Windows:
Run the following in R:

```
install.packages("devtools")
library(devtools)
devtools::install_github("protViz/fgcz.gsea.ora", build_vignettes = FALSE)
```
For installing package together with vignettes run:

```
devtools::install_github("protViz/fgcz.gsea.ora", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

# Vignettes:

https://rdrr.io/github/protViz/fgcz.gsea.ora/#vignettes




