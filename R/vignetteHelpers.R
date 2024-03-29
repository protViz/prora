.scriptCopyHelperVec <- function(runscripts, workdir = getwd(), packagename = "prora" ){
  res <- NULL
  for(scripts in runscripts){
    src_script <- file.path( find.package(packagename) , scripts )
    dest_script <- file.path(workdir ,basename(scripts))
    message("copy ", src_script, " to ", dest_script)

    if (!file.copy(src_script , dest_script, overwrite = TRUE)) {
      warning(paste("could not copy script file.", dest_script, sep = " "))
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

#' copy all files needed to generate a bfabric report
#'
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_bfabric_2grpGSEA <- function(workdir = getwd()){
  runscripts <- c("analyseSingle2GrpGSEA/VisualizeSingle.Rmd",
                  "analyseSingle2GrpGSEA/ErrorMessage.Rmd")
  .scriptCopyHelperVec(runscripts, workdir = workdir)
}

#' copy all files needed to generate a bfabric report for cluster profiling
#'
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_bfabric_ClusterProfiler <- function(workdir = getwd()){
  runscripts <- c("clusterProfiler/profileClusters_V2.Rmd")
  .scriptCopyHelperVec(runscripts, workdir = workdir)
}


#' copy multigroup analysis
#'
#' @param workdir directory where to copy file - default is current working directory.
#' @export
#'
copy_multigroup_GSEA <- function(workdir = getwd()){
  runscripts <- c("analyseMultiple2GrpOutputs/multigroupGSEA.Rmd",
                  "analyseMultiple2GrpOutputs/prepare.R")
  .scriptCopyHelperVec(runscripts, workdir = workdir)
}
