#' Calculate SpatialFeatures
#'
#' @description
#' This function takes a MoleculeExperiment (ME) object as input and calculates
#' new molecule-based spatial features according to the following segmentations:
#'
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Super-sector polygons
#' - Super-concentric polygons
#'
#' The function has 3 key steps, 1. adds new boundaries to the ME object
#' according to the subcellular and supercellular segmentations, 2. calculates
#' entropy for each of the cell and featureType combinations, and 3. combines
#' the entropy values into a SummarizedExperiment object.
#'
#' @param me A MoleculeExperiment (ME) object
#' @param featureTypes a character vector listing the types of featureTypes to include
#' @param k A numeric value indicating the number of polygons to calculate entropy (default 5)
#' @param nCores number of cores for parallel processing (default 1)
#' @param includeCounts logical (default FALSE) whether to include gene counts as features
#' @param ... arguments passed to loadBoundaries and EntropyMatrix
#' @return A SummarizedExperiment object containing a spatialFeatures assay and cell-level colData
#' @export
#' @examples
#' data(example_me)
#' se <- spatialFeatures(me)
#' se
spatialFeatures <- function(me,
                            featureTypes = c("sub_sector", "sub_concentric",
                                             "super_sector", "super_concentric"),
                            k = 5,
                            nCores = 1,
                            includeCounts = FALSE,
                            ...) {

  # step 1 load new boundaries
  me <- loadBoundaries(me, k = k, ...)

  # step 2 calculate entropies
  ent = EntropyMatrix(me, nCores = nCores, featureTypes = featureTypes, ...)

  # step 3 create SummarizedExperiment
  se = EntropySummarizedExperiment(ent, me, includeCounts = includeCounts, nCores = nCores)

  return(se)
}
