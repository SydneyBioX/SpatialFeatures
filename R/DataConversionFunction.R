#' Convert to Long Format
#' 
#' Converts a given entropy table to a long format.
#' 
#' @param me An object containing information about the experiment.
#' @param assayName Character string specifying the assay type.
#' 
#' @return A data frame in long format.
convert_to_long <- function(me, assayName) {
  mat = EntropyMatrix(me, assayName)
  df =  as.data.frame(mat) %>% 
    rownames_to_column("Gene_j") %>%
    gather(key = "Cell_i", value = "Value (Entropy)", -Gene_j) %>%
    mutate(Approach = assayName) %>%
    select("Value (Entropy)", Cell_i, Gene_j, Approach)
  return(df)
}