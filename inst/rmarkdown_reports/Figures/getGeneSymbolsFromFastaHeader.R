
getSymbolFromFasta <- function(fasta_header) {
  ref_protein_list <- data.frame(IDs = fasta_header[grepl("sp", fasta_header)]) %>%
    separate(col = IDs,
             sep = "_",
             into = c("begin", "end")) %>% 
    separate(col = begin,
             sep = "\\|",
             into = c("prefix","uniprotid","Symbol")) %>%
    select(Symbol) %>% 
    distinct
  return(ref_protein_list)
}
