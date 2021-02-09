
#' create url for data table
#' @export
DT_makeURLfor <- function(x, path){
  url <- paste0(path,x)
  html <- paste0("<a href='",url,"' target='_blank'>", x  ,"</a>")
  return(html)
}
