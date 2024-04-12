#' Load Boundaries to Molecule Experiment Object
#'
#' @description 
#' This function takes a Molecule Experiment (ME) object as input and enriches it 
#' with boundary data from a variety of assays. The currently supported assays are:
#'
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Sub-combo polygons
#' - Super-concentric polygons
#' - Super-combo polygons
#'
#' The function internally processes each assay, extracts respective boundary information, 
#' and updates the original ME object with this new data.
#'
#' @param me A Molecule Experiment (ME) object that needs its boundaries populated based 
#'        on the assays mentioned above.
#' @return An enriched Molecule Experiment (ME) object that includes boundary information 
#'         for the supported assays.
#' @export
#' @examples
#' data(example_me)
#' me_list <- loadBoundariesToMEforAll(example_me)

loadBoundariesToMEforAll <- function(me, k = 5) {
  
  regions_segments <- split_rectangle(me)
  
  assay_names <- c("sub_sector", "sub_concentric", 
                   "super_sector", "super_concentric")
  
  results_all <- SpatialFeatures::GenerateFeatureData(me, k = k)
  
  df_boundaries <- MoleculeExperiment::boundaries(me, assayName = "cell", flatten = TRUE)
  
  df_molecule <- MoleculeExperiment::molecules(me, assayName = "detected", flatten = TRUE)
  
  process_segment <- function(i) {
    results_list <- list()
    moleculesMEList_current <- list()
    boundariesMEList_current <- list()
    
    df_bdy = split_dataframe_by_region(df_boundaries, regions_segments)[[i]]
    
    if(nrow(df_bdy) == 0) {
      return("null")
    }
    
    boundariesMEList_current <- dataframeToMEList(df = df_bdy, dfType = "boundaries", assayName = "cell", 
                                                  sampleCol = "sample_id", factorCol = "segment_id", 
                                                  xCol = "x_location", yCol = "y_location", 
                                                  keepCols = "essential", scaleFactor = 1)
    
    for (assayName in assay_names) {
      df_tmp <- split_dataframe_by_region(results_all[[assayName]], regions_segments)[[i]]
      
      if (nrow(df_tmp) == 0) next
      
      df_tmp$segment_id <- df_tmp$area_id
      df_tmp$x_location <- df_tmp$x_section
      df_tmp$y_location <- df_tmp$y_section
      
      if (assayName == "super_sector") {
        coord_ranges <- df_tmp %>%
          dplyr::group_by(sample_id) %>%
          dplyr::summarise(
            xmin = min(x_location, na.rm = TRUE),
            xmax = max(x_location, na.rm = TRUE),
            ymin = min(y_location, na.rm = TRUE),
            ymax = max(y_location, na.rm = TRUE),
            .groups = "drop"
          )
        
        filtered_molecule <- df_molecule %>%
          dplyr::inner_join(coord_ranges, by = "sample_id") %>%
          dplyr::filter(x_location >= xmin & x_location <= xmax &
                          y_location >= ymin & y_location <= ymax)
        
        filtered_molecule$x_coords <- filtered_molecule$x_location
        filtered_molecule$y_coords <- filtered_molecule$y_location
        
        moleculesMEList_current <- dataframeToMEList(filtered_molecule, dfType = "molecules", assayName = "detected", 
                                                     sampleCol = "sample_id", factorCol = "features", xCol = "x_coords", 
                                                     yCol = "y_coords")
      }
      
      results_list[[assayName]] <- dataframeToMEList(df = df_tmp, dfType = "boundaries", assayName = assayName, 
                                                     sampleCol = "sample_id", factorCol = "segment_id", 
                                                     xCol = "x_location", yCol = "y_location", 
                                                     keepCols = "essential", scaleFactor = 1)
    }
    
    me_copy <- MoleculeExperiment(
      molecules = moleculesMEList_current,
      boundaries = boundariesMEList_current
    )
    
    for (assay in assay_names) {
      boundaries(me_copy, assay) <- results_list[[assay]]
    }
    
    return(me_copy)
  }
  
  me_list <- lapply(1:12, process_segment)
  
  return(me_list)
}


# loadBoundariesToME <- function(me, k = 8) {
#   
#   regions_segments <- split_rectangle(me)
#   
#   # List of assay names to process
#   assay_names <- c("sub_sector", "sub_concentric", "sub_combo", 
#                    "super_sector", "super_concentric", "super_combo")
#   
#   # Generate all the data using the GenerateFeatureData function
#   results_all <- GenerateFeatureData(me, k = k)
#   
#   # Read the molecule data
#   df_molecule <- MoleculeExperiment::molecules(me, assayName = "detected", flatten = TRUE)
#   
#   process_segment <- function(i) {
#     # List to store the results for each assay name
#     results_list <- list()
#     moleculesMEList_current <- list()
#     
#     for (assayName in assay_names) {
#       
#       # Split the dataframe based on region
#       df_tmp <- split_dataframe_by_region(results_all[[assayName]], regions_segments)[[i]]
#       
#       # Check if the data frame is empty
#       if (nrow(df_tmp) == 0) next
#       
#       # Modify column names based on the given logic
#       df_tmp$segment_id <- df_tmp$area_id
#       df_tmp$x_location <- df_tmp$x_section
#       df_tmp$y_location <- df_tmp$y_section
#       
#       # Filter the location of Molecules data
#       if (assayName == "super_sector") {
#         coord_ranges <- df_tmp %>%
#           dplyr::group_by(sample_id) %>%
#           dplyr::summarise(
#             xmin = min(x_location, na.rm = TRUE),
#             xmax = max(x_location, na.rm = TRUE),
#             ymin = min(y_location, na.rm = TRUE),
#             ymax = max(y_location, na.rm = TRUE),
#             .groups = "drop"
#           )
#         
#         filtered_molecule <- df_molecule %>%
#           dplyr::inner_join(coord_ranges, by = "sample_id") %>%
#           dplyr::filter(x_location >= xmin & x_location <= xmax &
#                           y_location >= ymin & y_location <= ymax)
#         
#         filtered_molecule$x_coords <- filtered_molecule$x_location
#         filtered_molecule$y_coords <- filtered_molecule$y_location
#         
#         moleculesMEList_current <- dataframeToMEList(filtered_molecule, dfType = "molecules", assayName = "detected", 
#                                                      sampleCol = "sample_id", factorCol = "features", xCol = "x_coords", 
#                                                      yCol = "y_coords")
#       }
#       
#       # Convert to ME List using dataframeToMEList function 
#       results_list[[assayName]] <- dataframeToMEList(df = df_tmp, dfType = "boundaries", assayName = assayName, 
#                                                      sampleCol = "sample_id", factorCol = "segment_id", 
#                                                      xCol = "x_location", yCol = "y_location", 
#                                                      keepCols = "essential", scaleFactor = 1)
#     }
#     
#     me_copy <- MoleculeExperiment(
#       molecules = moleculesMEList_current,
#       # boundaries = me@boundaries
#     )
#     
#     # Set the boundaries for each assay in me_copy
#     for (assay in assay_names) {
#       boundaries(me_copy, assay) <- results_list[[assay]]
#     }
#     
#     return(me_copy)
#   }
#   
#   # Use mclapply for parallel processing
#   me_copies <- mclapply(1:12, process_segment, mc.cores = 10)
#   
#   return(me_copies)
# }

# loadBoundariesToME <- function(me, k = 8) {
#   
#   regions_segments <- split_rectangle(me)
#   
#   # List of assay names to process
#   assay_names <- c("sub_sector", "sub_concentric", "sub_combo", 
#                    "super_sector", "super_concentric", "super_combo")
#   
#   # Generate all the data using the GenerateFeatureData function
#   results_all <- GenerateFeatureData(me, k = 8)
#   
#   me_copies <- list()
#   
#   # Dictionary to store the results for each assay name
#   results_list <- list()
#   
#   # Read the molecule data
#   df_molecule = MoleculeExperiment::molecules(me, assayName = "detected", flatten = TRUE)
#   
#   moleculesMEList <- list()
#   
#   # Create copies of me
#   for (i in 1:12) {
#     # Dictionary to store the results for each assay name
#     results_list <- list()
#     # Reset moleculesMEList for every iteration
#     moleculesMEList_current <- moleculesMEList
#     for (assayName in assay_names) {
#       
#       # Split the dataframe based on region
#       list_of_dfs <- split_dataframe_by_region(results_all[[assayName]], regions_segments)
#       df_tmp <- list_of_dfs[[i]]
#       
#       # Check if the data frame is empty
#       if (nrow(df_tmp) == 0) {
#         next  # skip to the next iteration if data frame is empty
#       }
#       
#       # Modify column names based on the given logic
#       df_tmp$segment_id <- df_tmp$area_id
#       df_tmp$x_location <- df_tmp$x_section
#       df_tmp$y_location <- df_tmp$y_section
#       
#       # Filter the location of Molucules data
#       if (assayName == "super_sector") {
#         # First, group by sample_id and then determine min/max coordinates
#         coord_ranges <- df_tmp %>%
#           group_by(sample_id) %>%
#           summarise(
#             xmin = min(x_location, na.rm = TRUE),
#             xmax = max(x_location, na.rm = TRUE),
#             ymin = min(y_location, na.rm = TRUE),
#             ymax = max(y_location, na.rm = TRUE),
#             .groups = "drop"
#           )
#         
#         # Filter the molecules based on these coordinates and sample_id
#         filtered_molecule <- df_molecule %>%
#           inner_join(coord_ranges, by = "sample_id") %>%
#           filter(x_location >= xmin & x_location <= xmax &
#                    y_location >= ymin & y_location <= ymax)
#         filtered_molecule$x_coords <- filtered_molecule$x_location
#         filtered_molecule$y_coords <- filtered_molecule$y_location
#         moleculesMEList_current <- dataframeToMEList(filtered_molecule, dfType = "molecules", assayName = "detected", 
#                                                      sampleCol = "sample_id", factorCol = "features", xCol = "x_coords", 
#                                                      yCol = "y_coords")
#       }
#       
#       # Convert to ME List using dataframeToMEList function 
#       results_list[[assayName]] <- dataframeToMEList(df = df_tmp, dfType = "boundaries", assayName = assayName, 
#                                                      sampleCol = "sample_id", factorCol = "segment_id", 
#                                                      xCol = "x_location", yCol = "y_location", 
#                                                      keepCols = "essential", scaleFactor = 1)
#     }
#     
#     me_copies[[i]] <- MoleculeExperiment(
#       molecules = moleculesMEList_current,
#       boundaries = me@boundaries
#     )
#     
#     # Load data into me using boundaries function
#     boundaries(me_copies[[i]], "sub_sector") <- results_list[["sub_sector"]]
#     boundaries(me_copies[[i]], "sub_concentric") <- results_list[["sub_concentric"]]
#     boundaries(me_copies[[i]], "sub_combo") <- results_list[["sub_combo"]]
#     boundaries(me_copies[[i]], "super_sector") <- results_list[["super_sector"]]
#     boundaries(me_copies[[i]], "super_concentric") <- results_list[["super_concentric"]]
#     boundaries(me_copies[[i]], "super_combo") <- results_list[["super_combo"]]
#   }
#   
#   return(me_copies)
# }

# loadBoundariesToME <- function(me, k = 8) {
#   
#   # List of assay names to process
#   assay_names <- c("sub_sector", "sub_concentric", "sub_combo", 
#                    "super_sector", "super_concentric", "super_combo")
#   
#   # Generate all the data using the GenerateFeatureData function
#   results_all <- GenerateFeatureData(me, k)
#   
#   # Dictionary to store the results for each assay name
#   results_list <- list()
#   
#   for (assayName in assay_names) {
#     
#     # Extract the appropriate dataset from results_all
#     df_tmp <- results_all[[assayName]]
#     
#     # Modify column names based on the given logic
#     df_tmp$segment_id <- df_tmp$area_id
#     df_tmp$x_location <- df_tmp$x_section
#     df_tmp$y_location <- df_tmp$y_section
#     
#     # Convert to ME List using dataframeToMEList function 
#     results_list[[assayName]] <- dataframeToMEList(df = df_tmp, dfType = "boundaries", assayName = assayName, 
#                                                    sampleCol = "sample_id", factorCol = "segment_id", 
#                                                    xCol = "x_location", yCol = "y_location", 
#                                                    keepCols = "essential", scaleFactor = 1)
#   }
#   
#   # Load data into me using boundaries function
#   boundaries(me, "sub_sector") <- results_list[["sub_sector"]]
#   boundaries(me, "sub_concentric") <- results_list[["sub_concentric"]]
#   boundaries(me, "sub_combo") <- results_list[["sub_combo"]]
#   boundaries(me, "super_sector") <- results_list[["super_sector"]]
#   boundaries(me, "super_concentric") <- results_list[["super_concentric"]]
#   boundaries(me, "super_combo") <- results_list[["super_combo"]]
#   
#   return(me)
# }