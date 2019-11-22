library(sigora)
library(SRMService)
library(tidyverse)
library(org.Hs.eg.db)

# packagedir <- path.package("SRMService")
# protein <- readr::read_tsv(file.path(packagedir, "samples/proteinGroups/proteinGroupsPullDown.txt"))

load("/Users/Lucas/Dropbox/FGCZ/SRMService/data/webGestaltExample.rda")

ref <-
  select(
    org.Hs.eg.db,
    columns = c("SYMBOL", "UNIPROT", "ENSEMBL"),
    keys = as.character(webGestaltExample$reference_list$rownames.quant_data.),
    keytype = "UNIPROT",
    multiVals = "first"
  )


targ <-
  select(
    org.Hs.eg.db,
    columns = c("SYMBOL"),
    keys = as.character(webGestaltExample$clusterIDs$colID[webGestaltExample$clusterIDs$clusterID ==
                                                             1]),
    keytype = "UNIPROT",
    multiVals = "first"
  )

colnames(ref) <- c("Uniprot", "Symbol")

## using precompiled GPS repository:
sigRes.ilreact <- sigora(queryList = targ$SYMBOL,
                         GPSrepo = kegH,
                         level = 3)


a1 <- genesFromRandomPathways(
  GPSrepo = kegH,
  np = 10,
  # From np random pathways
  ng = 150,
  # Draw nq genes at random
  seed = seed
)

reference_mapped_mapped <-
  dplyr::inner_join(data.frame(EntrezGene.ID = a1[[1]]), idmap)

# head(idmap)
# 
# reference_list <- getSymbolFromFasta(protein$`Protein IDs`)
# 
# reference_mapped <- inner_join(idmap, ref)
# 
# reference_mapped_mapped <- inner_join(reference_mapped, nciTable)

## user created GPS repository:
nciH <- makeGPS(pathwayTable = reference_mapped_mapped)
sigRes.ilnci <- sigora(queryList = targ$SYMBOL,
                       GPSrepo = nciH,
                       level = 3)

