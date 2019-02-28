# sigora wrapper functions

library(tidyverse)

getSymbolFromFasta <- function(df) {
  df %>%
    dplyr::filter(grepl(pattern = "sp", df$protein_Id)) %>%
    separate(col = protein_Id,
             sep = "_",
             into = c("begin", "end")) %>%
    separate(
      col = begin,
      sep = "\\|",
      into = c("prefix", "uniprotid", "Symbol")
    ) %>%
    dplyr::select(-prefix, -uniprotid, -end) %>%
    return()
}

sigoraWrappR <-
  function(input.file = "",
           fc_threshold = 0.5,
           fc_col = "",
           GPSrepos = kegH,
           df,
           db = "") {
    enriched <- df[df[, fc_col] >= fc_threshold,]
    sigora_res <-
      sigora(GPSrepo = GPSrepos,
             level = 5,
             queryList = enriched$Symbol)
    ora_res <-
      ora(geneList = enriched$Symbol,
          GPSrepo = GPSrepos)
    output <- list(
      sigora = sigora_res,
      ora = ora_res,
      fc_threshold = fc_threshold,
      GPSrepository = deparse(substitute(GPSrepos)),
      database = db,
      data = df[, c("Symbol", fc_col)],
      proteinsAfterFiltering = nrow(df[df[, fc_col] >= fc_threshold,])
    )
    return(output)
  }

# Function for generating background GPS repository for sigora and ora

makeGPS_wrappR <- function(ids, target = "KEGG") {
  if (target == "KEGG") {
    gp_db <- org.Hs.eg.db
    target_column <- "PATH"
    pn_table <- sigora::kegH$pathwaydescriptions
    colnames(pn_table) <- c("pathwayID", "pathwayName")
    pn_table$pathwayID <-
      substr(pn_table$pathwayID, start = 4, stop = 8)
    gp_table <- AnnotationDbi::mapIds(
      gp_db,
      keys = ids,
      keytype = "SYMBOL",
      column = target_column,
      multiVals = "CharacterList"
    ) %>% unlist %>% data.frame(Symbol = names(.), pathwayID = .)
  } else {
    if (target == "reactome") {
      gp_db <- reactome.db
      pn_db <- reactome.db
      ids <- AnnotationDbi::mapIds(
        # First map to ENTREZID
        org.Hs.eg.db,
        keys = ids,
        keytype = "SYMBOL",
        column = "ENTREZID",
        multiVals = "CharacterList"
      ) %>% unlist %>% na.omit
      target_column <- "PATHID"
      k0 = "ENTREZID"
      k1 = "PATHID"
      k2 = "PATHNAME"
    } else
      if (target == "GO") {
        gp_db <- org.Hs.eg.db
        pn_db <- GO.db
        target_column = "GO"
        k0 = "SYMBOL"
        k1 = "GOID"
        k2 = "TERM"
      } else
        return(message("Specify a valid target database (KEGG, reactome, GO)"))

    gp_table <- mapIds(
      gp_db,
      keys = ids,
      keytype = k0,
      column = target_column,
      multiVals = "CharacterList"
    ) %>% unlist %>% data.frame(Symbol = names(.), pathwayID = .)

    pn_table <- mapIds(
      pn_db,
      keys = as.character(gp_table$pathwayID),
      multiVals = "CharacterList",
      keytype = k1,
      column = k2
    ) %>% unlist %>% data.frame(pathwayID = names(.), pathwayName = .)
  }

  mkGPStable <- inner_join(pn_table, gp_table) %>% distinct

  colnames(mkGPStable) <- c("pathwayId", "pathwayName", "gene")

  out <- makeGPS(mkGPStable)

  return(out)
}
