#' Load Boundaries to Molecule Experiment Object
#'
#' @description
#' This function takes a Molecule Experiment (ME) object as input and enriches
#' it with boundary data from a variety of assays. The currently supported
#' assays are:
#'
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Super-sector polygons
#' - Super-concentric polygons
#'
#' The function internally processes each assay, extracts respective boundary
#' information, and updates the original ME object with this new data.
#'
#' @param me A MoleculeExperiment (ME) object
#' @param ... arguments passed to GenerateFeatureData
#' @return An enriched Molecule Experiment (ME) object that includes boundary
#' information for the supported assays.
#' @importFrom MoleculeExperiment dataframeToMEList
#' @importFrom MoleculeExperiment boundaries
#' @export
#' @examples
#' data(example_me)
#' me <- loadBoundaries(me)
#' me
loadBoundaries <- function(me,...) {

  featureData <- GenerateFeatureData(me, ...)

  featureData <- lapply(featureData,
                        function(x) {x$segment_id <- x$area_id; return(x)})

  bds_list <- lapply(featureData, MoleculeExperiment::dataframeToMEList,
                    dfType = "boundaries", assayName = "cell",
                    sampleCol = "sample_id", factorCol = "area_id",
                    xCol = "x_section", yCol = "y_section",
                    keepCols = "essential", scaleFactor = 1
  )

  for (assay in names(featureData)) {
    MoleculeExperiment::boundaries(me, assay) <- bds_list[[assay]]
  }

  return(me)
}
