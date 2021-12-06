
#' create url for data table
#' @param x eg. uniprot id's P94485
#' @param path e.g. https://www.uniprot.org/uniprot/
#' @export
DT_makeURLfor <- function(x, path){
  url <- paste0(path,x)
  html <- paste0("<a href='",url,"' target='_blank'>", x  ,"</a>")
  return(html)
}
