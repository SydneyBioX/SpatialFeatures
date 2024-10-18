#' Create a SingleCellExperiment Object from a List of Data Frames
#'
#' This function takes a list of data frames and a Molecule Experiment object
#' and constructs a \linkS4class{SingleCellExperiment} object with the
#' specified data structure.
#'
#' @param df_list A list of data frames, each containing assay data.
#' @param me A \linkS4class{MoleculeExperiment} object.
#' @param includeCounts logical (default FALSE) whether to include gene
#' counts as features
#' @param concatenateFeatures logical whether to concatenate all the features
#' into a single assay (default FALSE). If FALSE the output SCE object has
#' multiple assays
#' @param nCores Number of cores (default 1)
#'
#' @return A SingleCellExperiment object.
#'
#' @export
#' @import SingleCellExperiment
#' @import SpatialFeatures
#' @importFrom MoleculeExperiment countMolecules
#' @examples
#' data(example_me)
#' me <- loadBoundaries(me)
#' ent <- EntropyMatrix(me, c("subsector", "subconcentric",
#' "supersector", "superconcentric"), nCores = 1)
#' se <- EntropySingleCellExperiment(ent, me)
#' se
EntropySingleCellExperiment <- function(df_list,
                                        me,
                                        includeCounts = FALSE,
                                        concatenateFeatures = FALSE,
                                        nCores = 1) {

  # 1. Assay Data: Using countMolecules function to get assay data.
  # Creating the assay_data
  assay_data <- make_assay_data(df_list,
                                concatenateFeatures = concatenateFeatures)

  if (includeCounts) {
    # Generating the genecount
    genecount <- MoleculeExperiment::countMolecules(me,
                                                    moleculesAssay = "detected",
                                                    boundariesAssay = "cell",
                                                    matrixOnly = TRUE,
                                                    nCores = nCores)

    if (concatenateFeatures) {

      # Adding a prefix to the row names of genecount to differentiate it
      rownames(genecount) <- paste("genecount", rownames(genecount), sep = "_")

      genecount <- as.matrix(genecount)

      # Rbinding the assay_data and genecount together
      assay_data <- rbind(assay_data, genecount)

    } else {

      assay_data[["genecount"]] <- genecount

    }
  }

  if (concatenateFeatures) {

    # 2. Row Data
    rowData <- data.frame(
      FeatureCategory = gsub("_.*", "", rownames(assay_data)),
      FeatureGene = gsub(".*_", "", rownames(assay_data))
    )

  } else {
    rowData <- NULL
  }

  # 3. Column Data
  cell_df <- extract_boundaries_and_centroids(me)[[2]]
  colnames(cell_df) <- c("x_location", "y_location",
                         "Cell", "Sample_id", "x_central", "y_central")
  cell_df_list <- split(cell_df[, c('x_location', 'y_location')], cell_df$Cell)
  colData <- unique(cell_df[, c("Cell", "Sample_id", "x_central", "y_central")])
  colData$boundaries <- I(cell_df_list)
  colData$Cell <- paste0(colData$Sample_id, ".", colData$Cell)

  if (concatenateFeatures) {
    se <- SingleCellExperiment::SingleCellExperiment(
      assays = list(spatialFeatures = as.matrix(assay_data)),
      rowData = rowData,
      colData = colData
    )
  } else {

    rnames <- gsub(".*_", "", rownames(assay_data[[1]]))
    cnames <- gsub(".*_", "", colnames(assay_data[[1]]))

    assay_data <- lapply(assay_data, function(x){
      rownames(x) <- rnames
      colnames(x) <- cnames
      return(x)
    })

    se <- SingleCellExperiment::SingleCellExperiment(
      assays = assay_data,
      rowData = rowData,
      colData = colData
    )
  }

  return(se)
}
