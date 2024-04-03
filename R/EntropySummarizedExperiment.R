#' Create a SummarizedExperiment Object from a List of Data Frames
#'
#' This function takes a list of data frames and a Molecule Experiment object and
#' constructs a SummarizedExperiment object with the specified data structure.
#'
#' @param df_list A list of data frames, each containing assay data.
#' @param me A Molecule Experiment object.
#'
#' @return A SummarizedExperiment object.
#'
#' @examples
#' \dontrun{
#' # Assuming `em_list_full` is your list of data frames and `me` is your Molecule Experiment object
#' se <- createSummarizedExperiment(df_list = em_list_full, me = me)
#' }
#' @export
#' @import SummarizedExperiment
#' @import SpatialFeatures
EntropySummarizedExperiment <- function(df_list, me) {
  
  # 1. Assay Data: Using countMolecules function to get assay data.
  # Creating the assay_data
  assay_data <- do.call(rbind, lapply(names(df_list), function(assayName) {
    df <- df_list[[assayName]]
    rownames(df) <- paste(assayName, rownames(df), sep="_")
    return(df)
  }))
  
  # Generating the genecount
  genecount <- as.data.frame(as.matrix(MoleculeExperiment::countMolecules(me, 
                                                                          moleculesAssay = "detected", 
                                                                          boundariesAssay = "cell", 
                                                                          matrixOnly = TRUE, 
                                                                          nCores = 10)))
  
  # Adding a prefix to the row names of genecount to differentiate it
  rownames(genecount) <- paste("genecount", rownames(genecount), sep="_")
  
  # Rbinding the assay_data and genecount together
  assay_data <- rbind(assay_data, genecount)
  
  # 2. Row Data
  rowData <- data.frame(
    FeatureCategory = gsub("_.*", "", rownames(assay_data)),
    FeatureGene = gsub(".*_", "", rownames(assay_data))
  )
  
  # Translate "sub" and "super" prefixes to "Subcellular" and "Supercellular" respectively
  rowData$FeatureCategory <- gsub("^sub", "Subcellular", rowData$FeatureCategory)
  rowData$FeatureCategory <- gsub("^super", "Supercellular", rowData$FeatureCategory)
  rowData$FeatureCategory <- gsub("^genecount", "Genecount", rowData$FeatureCategory)
  
  # 3. Column Data
  cell_df <- SpatialFeatures::extract_boundaries_and_centroids(me)[[2]]
  # Renaming the columns for clarity and ease of use
  colnames(cell_df) <- c("x_location", "y_location", "Cell", "Sample_id", "x_central", "y_central")
  # Now, create numericList for each cell
  cell_df_list <- split(cell_df[, c('x_location', 'y_location')], cell_df$Cell)
  # Extract unique rows
  colData <- unique(cell_df[, c("Cell", "Sample_id", "x_central", "y_central")])
  # Combine x_central, y_central, and numericList
  colData$boundaries <- I(cell_df_list)
  # Convert Cell column to a format compatible with the columns in assay_data
  colData$Cell <- paste0(colData$Sample_id, ".", colData$Cell)
  
  # Create the SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(assay_data)),
    rowData = rowData,
    colData = colData
  )
  
  return(se)
}
