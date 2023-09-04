#' Bind Rows of Multiple Methods into a Single Dataframe
#'
#' This function binds rows of multiple processing methods into a single 
#' dataframe in long format.
#' 
#' @param me A data object that contains the input data.
#' @param methods A character vector containing the names of the methods to be 
#'   processed. They must be among "sub-sector", "sub-concentric", "sub-combo", 
#'   "super-concentric", and "super-combo".
#'
#' @return A dataframe in long format.
#' @export
BindRows <- function(me, methods) {
  # Check for valid methods
  valid_methods <- c("sub-sector", "sub-concentric", 'sub-combo', "super-concentric", "super-combo")
  if(any(!methods %in% valid_methods)) {
    stop("One or more invalid methods provided!")
  }
  
  # Convert each method to long format and bind
  dfs <- lapply(methods, function(method) {
    convert_to_long(me, method)
  })
  
  # Bind rows
  long_df <- do.call(rbind, dfs)
  return(long_df)
}