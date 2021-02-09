

remotes::install_bioc("AnnotationHub")
library(AnnotationHub)
ah <- AnnotationHub::AnnotationHub()
orgdb <- AnnotationHub::query(ah, c("OrgDb", "maintainer@bioconductor.org"))

specODB <- orgdb[[grep("Homo",orgdb$species)]]
egid <- AnnotationDbi::keys(specODB, "ENTREZID")
select(specODB, egid, c("SYMBOL", "GENENAME", "UNIPROT"), "ENTREZID") %>% dim()

