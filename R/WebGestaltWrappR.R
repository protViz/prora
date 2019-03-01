#' Wrapper for WebGestaltR overrepresentation analysis in SRMService
#' @export webGestaltWrapper webGestaltWrapper
#' @param WebGestaltWrappR WebGestaltWrappR
#' @import dplyr quantable WebGestaltR SRMService
#' @return list with results
#' @param grp2 grp2 S4 object
#' @param enrichDatabase Database used for over representation analysis
#' @param organism organism proteins in the grp2 object derive from
#' @param nrNas number of allowed NAs per protein. Protein matrix will be filtered using this threshold.
#' @param se_threshold standard error of protein expressions over all conditions. Will be used as a filtering threshold.

webGestaltWrapper <- function(grp2, enrichDatabase, organism, se_threshold, work_dir = "output", nrNas, method) {

  # Additional functions ------------------------------------------------------

  write_interesting_geneFile <- function(ID, df, output.dir) {
    df %>%
      dplyr::filter(clusterID == ID) %>%
      dplyr::select(colID) %>%
      write_tsv(
        path = paste(output.dir, "/protein_cluster_", ID, ".txt", sep = ""),
        col_names = F
      )
    message(paste("protein_cluster_",
                  ID,
                  ".txt was written to ",
                  output.dir,
                  sep = ""))
  }

  aggregate_results <- function(output) {
    makedataframe <- function(ll){ if (is.null(ll$enrichResult) || grepl("ERROR", ll$enrichResult))
      return(NULL)
      else {
        tmp <- as.data.frame(ll$enrichResult)
        tmp %>%
          dplyr::mutate(file.origin = parse_number(tools::file_path_sans_ext(basename(ll$filename))))
      }
    }

    res <- lapply(output, makedataframe) %>%
      do.call(rbind, .)
    return(res)
  }

  preprocess_df <- function(grp2) {
    quant_data <- grp2$getNormalized()$data
    quant_data <- quant_data[grepl("sp", row.names(quant_data)), ]
    ref_protein_list <- data.frame(IDs = row.names(quant_data)) %>%
      tidyr::separate(col = IDs,
               sep = "\\|",
               into = c("prefix","uniprotid","proteinname")) %>%
      dplyr::select(uniprotid)
    row.names(quant_data) <- make.unique(ref_protein_list$uniprotid)
    return(quant_data)
  }

  se_filter <- function(quant_data, se_threshold, nrNas) {
    se_observed <- apply(quant_data, 1, sd, na.rm = TRUE)
    quant_data_filtered <- quant_data[se_observed > se_threshold, ]
    idx <- which(apply(quant_data_filtered, 1, function(x) length(which(is.na(x))))<nrNas)
    quant_data <- quant_data_filtered[idx,]
    return(quant_data_filtered)
  }

  getProteinClustering <- function(quant_data_filtered, method, nclust) {
    clustering <-
      quantable::simpleheatmap3(
        scale(t(quant_data_filtered), scale = FALSE),
        labCol = row.names(quant_data_filtered),
        plot = FALSE,
        nrOfClustersCol = nclust,
        method = method
      )$Col
    return(clustering)
  }

  quant_data <- preprocess_df(grp2)


  # Reference list ----------------------------------------------------------

  reference_list <- data.frame(rownames(quant_data))
  write_tsv(reference_list, file.path(work_dir, "referencelist.txt"), col_names = F)


  # Protein cluster lists ---------------------------------------------------

  quant_data_filtered <- se_filter(quant_data, se_threshold, nrNas)

  clustering <- getProteinClustering(quant_data_filtered, method, grp2$getNumberOfClusters())

  sapply(
    unique(clustering$clusterID),
    write_interesting_geneFile,
    df = clustering,
    output.dir = "output/ORA_inputFiles"
  )


  # WebGestaltR call --------------------------------------------------------

  output <-
    WebGestaltR::WebGestaltRBatch(
      enrichMethod = "ORA",
      organism = organism,
      enrichDatabase = enrichDatabase,
      interestGeneFolder = "output/ORA_inputFiles",
      referenceGeneFile = "output/referencelist.txt",
      is.output = FALSE,
      interestGeneType = "uniprotswissprot",
      referenceGeneType =  "uniprotswissprot",
      outputDirectory = "output/"
    )

  aggregated_results <- aggregate_results(output)

    # Output list -------------------------------------------------------------

    config <- list(
      normalisedData_unfiltered = quant_data,
      normalisedData_filtered = quant_data_filtered,
      numberOfProteinClusters = grp2$getNumberOfClusters(),
      enrichDatabase = enrichDatabase,
      organism = organism,
      webgestaltResults = aggregated_results,
      se_threshold = se_threshold,
      webgestaltList = output,
      clusterIDs = clustering,
      nrNas = nrNas,
      method = method,
      reference_list = reference_list
    )
    return(config)
}


