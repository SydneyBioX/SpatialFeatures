#' make assay data
#'
#' @param df_list A list of data frames
#' @param concatenateFeatures logical whether to concatenate
#' features (default FALSE)
#' @return if concatenateFeatures == TRUE, A single data frame,
#' otherwise a list containing data frames
make_assay_data = function(df_list, concatenateFeatures = FALSE) {

  assay_data_list <- lapply(names(df_list), function(assayName) {
    df <- df_list[[assayName]]
    rownames(df) <- paste(assayName, rownames(df), sep = "_")
    return(df)
  })
  names(assay_data_list) <- names(df_list)

  if (concatenateFeatures) {
    return(do.call(rbind, assay_data_list))
  } else {
    return(assay_data_list)
  }

}
