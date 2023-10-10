#' Convert Entropy Matrix to Long Format Dataframe
#'
#' @description 
#' This function converts an entropy matrix from a Molecule Experiment object to a long format dataframe.
#'
#' @param me A Molecule Experiment object.
#' @param assayName A character string specifying the assay type. Supported values include 
#' "sub-sector", "sub-concentric", "sub-combo", "super_sector", "super-concentric", and "super-combo".
#'
#' @return A dataframe in long format with columns "Value", "Cell", "Gene", and "Approach".
#' @export
#' @examples
#' \dontrun{
#' # Assuming `data_obj` is your Molecule Experiment object
#' long_df = convert_to_long(data_obj, assayName = "sub_sector")
#' }
ConverttoLong <- function(me, assayName, ...) {
  mat <- EntropyMatrix(me, assayName, ...)
  df <- as.data.frame(mat) %>% 
    rownames_to_column("Gene") %>%
    gather(key = "Cell", value = "Value", -Gene) %>%
    mutate(Approach = assayName) %>%
    select("Value", Cell, Gene, Approach)
  return(df)
}