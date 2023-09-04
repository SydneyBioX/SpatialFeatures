#' Load Boundaries to ME Data Object
#'
#' This function loads boundaries for multiple assay names into the specified 
#' ME data object.
#' 
#' @param me A data object to which the boundaries will be added.
#'
#' @return The updated ME data object with added boundaries.
#' @export
loadBoundariesToME <- function(me) {
  
  # List of assay names to process
  assay_names <- c("sub-sector", "sub-concentric", "sub-combo", 
                   "super-concentric", "super-combo")
  
  # Dictionary to store the results for each assay name
  results_list <- list()
  
  for (assayName in assay_names) {
    
    # Generate the data using the GenerateFeatureData function
    df_tmp <- GenerateFeatureData(me, assayName)
    
    # Modify column names based on the given logic (need to update after debug of Function dataframeToMEList)
    df_tmp$segment_id <- df_tmp$area_id
    df_tmp$x_location <- df_tmp$x_section
    df_tmp$y_location <- df_tmp$y_section
    
    # Convert to ME List using dataframeToMEList function 
    results_list[[assayName]] <- dataframeToMEList(df = df_tmp, dfType = "boundaries", assayName = assayName, 
                                                   sampleCol = "sample_id", factorCol = "segment_id", 
                                                   xCol = "x_location", yCol = "y_location", 
                                                   keepCols = "essential", scaleFactor = 1)
  }
  
  # Load data into me using boundaries function
  boundaries(me, "sub-sector") <- results_list[["sub-sector"]]
  boundaries(me, "sub-concentric") <- results_list[["sub-concentric"]]
  boundaries(me, "sub-combo") <- results_list[["sub-combo"]]
  boundaries(me, "super-concentric") <- results_list[["super-concentric"]]
  boundaries(me, "super-combo") <- results_list[["super-combo"]]
  
  return(me)
}