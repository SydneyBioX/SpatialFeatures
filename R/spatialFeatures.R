#' Calculate SpatialFeatures from a MoleculeExperiment object
#'
#' @description
#' This function takes a \linkS4class{MoleculeExperiment} (ME) object as input
#' and calculates new molecule-based spatial features according to the
#' following segmentations:
#'
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Super-sector polygons
#' - Super-concentric polygons
#'
#' The function has 3 key steps: 1. adds new boundaries to the ME object
#' according to the subcellular and supercellular segmentations, 2. calculates
#' entropy for each of the cell and featureType combinations, and 3. combines
#' the entropy values into a \linkS4class{SingleCellExperiment} object.
#'
#' @param me A MoleculeExperiment (ME) object
#' @param featureTypes A character string specifying the feature type.
#' Supported values include "subsector", "subconcentric", "supersector", and
#' "superconcentric".
#' @param k A numeric value indicating the number of polygons to calculate
#' entropy (default 5)
#' @param nCores number of cores for parallel processing (default 1)
#' @param includeCounts logical (default FALSE) whether to include gene
#' counts as features
#' @param concatenateFeatures logical whether to concatenate all the features
#' into a single assay (default FALSE). If FALSE the output SE object has
#' multiple assays
#' @param ... arguments passed to \code{\link{loadBoundaries}} and
#' \code{\link{EntropyMatrix}}
#' @return A SingleCellExperiment object containing a spatialFeatures assay
#' and cell-level colData
#' @export
#' @examples
#' data(example_me)
#' se <- spatialFeatures(me)
#' se
spatialFeatures <- function(me,
                            featureTypes = c("subsector", "subconcentric",
                                             "supersector", "superconcentric"),
                            k = 5,
                            nCores = 1,
                            includeCounts = FALSE,
                            concatenateFeatures = FALSE,
                            ...) {

  # step 0 check me is valid, contains boundaries assays
  check_if_me_boundaries(me)

  # step 1 load new boundaries
  me <- loadBoundaries(me, k = k, ...)

  # step 2 calculate entropies
  ent <- EntropyMatrix(me, nCores = nCores, featureTypes = featureTypes, ...)

  # step 3 create SingleCellExperiment
  se <- EntropySingleCellExperiment(ent, me, includeCounts = includeCounts,
                                    concatenateFeatures = concatenateFeatures,
                                    nCores = nCores)

  return(se)
}
