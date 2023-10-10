#' Bind Rows of Converted Entropy Matrix Data
#'
#' @description 
#' This function binds rows of entropy matrix data from a Molecule Experiment object for multiple methods into a single dataframe.
#'
#' @param me A Molecule Experiment object.
#' @param methods A character vector specifying the assay types to process. Supported values include 
#' "sub-sector", "sub-concentric", "sub-combo", "super-concentric", and "super-combo".
#'
#' @return A dataframe with bound rows from the provided methods.
#' @export
#' @examples
#' \dontrun{
#' # Assuming `data_obj` is your Molecule Experiment object
#' bound_df = BindRows(data_obj, methods = c("sub-sector", "sub-concentric"))
#' }
BindRows <- function(me, methods, ...) {

  # Check for valid methods
  valid_methods <- c("sub_sector", "sub_concentric", 'sub_combo', "super_sector","super_concentric", "super_combo")
  if (any(!methods %in% valid_methods)) {
    stop("One or more invalid methods provided!")
  }

  # Convert each method to long format and bind
  dfs <-mclapply(methods, function(method) {
    SpatialFeatures::ConverttoLong(me, method, ...)
  }, mc.cores = 10)

  # Bind rows
  long_df <- do.call(rbind, dfs)

  return(long_df)
}