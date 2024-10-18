#' Check that a MoleculeExperiment object is both valid and contains boundaries
#'
#' @param me A MoleculeExperiment object
#'
#' @importFrom MoleculeExperiment showBoundaries
#' @return Will return an error if me is not valid and containing boundaries
check_if_me_boundaries <- function(me) {

  MoleculeExperiment:::.check_if_me(me)

  if (is.null(me@boundaries)) {
    stop("MoleculeExperiment object must have at least one boundaries assay")
  }

}
