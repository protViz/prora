library(tidyverse)
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

x <- fgczgseaora::exampleUniprotIDs # 2292 unique uniprot ids from proteomics experiment

x <- fgczgseaora::exampleContrastData
head(x)

x <- x %>% dplyr::select(protein_Id) %>% separate(protein_Id, c("sp","uniID","genesymbol")) %>% pull(uniID)
length(unique(x))
quantable::write.vector(unique(x),file="forUniprot.txt")

Symbol_mapped <- mapIds( # substantial loss when mapping from symbol to GO
  x = org.Hs.eg.db,
  keys = x,
  column = "SYMBOL",
  keytype = "UNIPROT",
  multiVals = "CharacterList"
) %>% unlist() %>% data.frame(Symbol = names(.), ID = .) %>% distinct(Symbol, .keep_all = TRUE) %>% na.omit()

head(Symbol_mapped)
(1-nrow(Symbol_mapped)/length(x))*100 # percent loss in IDs
#2,288 out of 2,304 identifiers from UniProtKB AC/ID were successfully mapped to 2,307 Gene name IDs.


GO_mapped <- mapIds( # substantial loss when mapping from symbol to GO
  x = org.Hs.eg.db,
  keys = x,
  column = "GO",
  keytype = "UNIPROT",
  multiVals = "CharacterList"
) %>% unlist() %>% data.frame(Symbol = names(.), ID = .) %>% distinct(Symbol, .keep_all = TRUE) %>% na.omit()
(1-nrow(GO_mapped)/length(x))*100 # percent loss in IDs

entrezid_mapped <- mapIds( # substantial loss when mapping from symbol to ENTREZID
  x = org.Hs.eg.db,
  keys = x,
  column = "ENTREZID",
  keytype = "UNIPROT",
  multiVals = "CharacterList"
) %>% unlist() %>% data.frame(Symbol = names(.), ID = .) %>% distinct(Symbol, .keep_all = TRUE) %>% na.omit()
(1-nrow(entrezid_mapped)/length(x))*100 # percent loss in IDs
#2,273 out of 2,304 identifiers from UniProtKB AC/ID were successfully mapped to 2,305 Entrez Gene (GeneID) IDs.
(1- (2273/ 2304))*100

path_mapped <- mapIds( # substantial loss of IDs when mapping to KEGG
  x = org.Hs.eg.db,
  keys = x,
  column = "PATH",
  keytype = "UNIPROT",
  multiVals = "CharacterList"
) %>% unlist() %>% data.frame(Symbol = names(.), ID = .) %>% distinct(Symbol, .keep_all = TRUE) %>% na.omit()
(1-nrow(path_mapped)/length(x))*100 # percent loss in IDs
# 2,268 out of 2,304 identifiers from UniProtKB AC/ID were successfully mapped to 2,298 KEGG IDs.

reactome_mapped <- mapIds( # again a loss of IDs when mapping from entrez to reactomeid
  x = reactome.db,
  keys = as.character(entrezid_mapped$ID),
  column = "REACTOMEID",
  keytype = "ENTREZID",
  multiVals = "CharacterList"
) %>% unlist() %>% data.frame(Symbol = names(.), ID = .) %>% distinct(Symbol, .keep_all = TRUE) %>% na.omit()

dim(reactome_mapped)
(1-nrow(reactome_mapped)/nrow(entrezid_mapped))*100 # percent loss in IDs from second to third mapping
(1-nrow(reactome_mapped)/length(x))*100 # percent loss in IDs total for reactome

quantable::write.vector(fgczgseaora::exampleSymbols, file="forUniprot.txt")
