#' make assay data
#'
#' @param df_list A list of data frames
#' @return A single data frame
make_assay_data = function(df_list) {
  do.call(rbind, lapply(names(df_list), function(assayName) {
    df <- df_list[[assayName]]
    rownames(df) <- paste(assayName, rownames(df), sep="_")
    return(df)
  }))
}
