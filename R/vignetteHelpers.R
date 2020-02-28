.scriptCopyHelperVec <- function(runscripts, workdir = getwd(), packagename = "fgcz.gsea.ora" ){
  res <- NULL
  for(scripts in runscripts){
    src_script <- file.path( find.package(packagename) , scripts )
    dest_script <- file.path(workdir ,basename(scripts))
    message("copy ", src_script, " to ", dest_script)

    if(!file.copy(src_script , dest_script, overwrite = TRUE)){
      warning(paste("could not copy script file.", dest_script, sep=" "))
    }else{
      res <- c(res, dest_script )
    }
  }
  message(paste("your working directory now should contain: ", length(res) , "new files :\n",sep=" "))
  return(res)
}


#' copy all files needed to generate GSEA report.
#'
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_gsea_report <- function(workdir = getwd()){

  runscripts <- c("fgcz_formatting/fgcz_header.html",
                  "fgcz_formatting/fgcz_footer.html",
                  "fgcz_formatting/fgcz.css",
                  "fgcz_formatting/fgcz_banner.png",
                  "rmarkdown_reports/GSEA_Results_Overview.Rmd")

  .scriptCopyHelperVec(runscripts, workdir = workdir)
}
#' copy all files needed to generate ORA report
#'
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_ora_report <- function(workdir = getwd()){

  runscripts <- c("fgcz_formatting/fgcz_header.html",
                  "fgcz_formatting/fgcz_footer.html",
                  "fgcz_formatting/fgcz.css",
                  "fgcz_formatting/fgcz_banner.png",
                  "rmarkdown_reports/ORA_Results_Overview.Rmd")

  .scriptCopyHelperVec(runscripts, workdir = workdir)
}
