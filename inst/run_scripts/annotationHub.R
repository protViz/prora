ah <- AnnotationHub()
orgdb <- query(ah, c("OrgDb", "maintainer@bioconductor.org"))
class(orgdb)
length(orgdb)
names(orgdb)
length(orgdb$species)
orgdb1 <- query(ah, c("OrgDb", "maintainer@bioconductor.org"))[[1]]
class(orgdb1)

orgdb1
egid <- keys(orgdb1, "ENTREZID")
select(orgdb1, egid, c("SYMBOL", "GENENAME", "UNIPROT"), "ENTREZID")
