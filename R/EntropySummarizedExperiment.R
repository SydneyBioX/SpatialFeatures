#' Create a SummarizedExperiment Object from a List of Data Frames
#'
#' This function takes a list of data frames and a Molecule Experiment object and
#' constructs a SummarizedExperiment object with the specified data structure.
#'
#' @param df_list A list of data frames, each containing assay data.
#' @param me A Molecule Experiment object.
#' @param includeCounts logical (default FALSE) whether to include gene counts as features
#' @param nCores Number of cores (default 1)
#'
#' @return A SummarizedExperiment object.
#'
#' @export
#' @import SummarizedExperiment
#' @import SpatialFeatures
#' @importFrom MoleculeExperiment countMolecules
#' @examples
#' data(example_me)
#' me <- loadBoundaries(me)
#' ent <- EntropyMatrix(me, c("sub_sector", "sub_concentric", "super_sector", "super_concentric"), nCores = 1)
#' se <- EntropySummarizedExperiment(ent, me)
#' se
EntropySummarizedExperiment <- function(df_list, me, includeCounts = FALSE, nCores = 1) {

  # 1. Assay Data: Using countMolecules function to get assay data.
  # Creating the assay_data
  assay_data <- make_assay_data(df_list)

  if (includeCounts) {
  # Generating the genecount
  genecount <- as.data.frame(as.matrix(MoleculeExperiment::countMolecules(me,
                                                                          moleculesAssay = "detected",
                                                                          boundariesAssay = "cell",
                                                                          matrixOnly = TRUE,
                                                                          nCores = nCores)))

  # Adding a prefix to the row names of genecount to differentiate it
  rownames(genecount) <- paste("genecount", rownames(genecount), sep="_")

  # Rbinding the assay_data and genecount together
  assay_data <- rbind(assay_data, genecount)
  }

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
  cell_df <- extract_boundaries_and_centroids(me)[[2]]
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
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(spatialFeatures = as.matrix(assay_data)),
    rowData = rowData,
    colData = colData
  )

  return(se)
}
