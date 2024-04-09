#' Generate Feature Data Based on Assay Type Combining Various Polygon Assays
#'
#' @description 
#' This comprehensive function processes the given data (`me`) and returns feature 
#' data based on multiple assay types. The supported assays encapsulated within 
#' this function include:
#' 
#' - Sub-sector polygons
#' - Sub-concentric polygons
#' - Sub-combo polygons
#' - Super-concentric polygons
#' - Super-combo polygons
#'
#' These assays represent various spatial subdivisions or transformations, each 
#' addressing unique feature characteristics. The results from each assay are 
#' combined and returned as a list of dataframes.
#'
#' @param me The input data object. This object contains spatial information, 
#'        typically regarding segments or regions of interest.
#' @param k A numeric value indicating the scaling factor for concentric polygons. 
#'        Defaults to 8.
#' @return A list containing dataframes for each assay type:
#' - `sub_sector`: Feature data for sub-sector polygons.
#' - `sub_concentric`: Feature data for sub-concentric polygons.
#' - `sub_combo`: Feature data for sub-combo polygons.
#' - `super_concentric`: Feature data for super-concentric polygons.
#' - `super_combo`: Feature data for super-combo polygons.
#' @export
#' @importFrom dplyr group_by mutate arrange select ungroup rowwise
#' @importFrom purrr map_dfr
#' @importFrom parallel mclapply
#' @examples
#' data(example_me)
#' feature_data_list <- GenerateFeatureData(example_me)
#' sub_sector_data <- feature_data_list$sub_sector

GenerateFeatureData <- function(me, k = 5) {
  results <- extract_boundaries_and_centroids(me)
  df_circle = results$df_circle
  
  # For sub-sector
  sub_sector <- df_circle %>%
    distinct(segment_id, x_location, y_location, .keep_all = TRUE) %>%
    mutate(angle = atan2(y_location - y_central, x_location - x_central) + pi) %>%
    split(.$segment_id) %>%
    lapply(create_sectors) %>%
    bind_rows() %>%
    dplyr::rename(x_section = x_location_sector, y_section = y_location_sector, area_id = sector_id) 
  
  # Common for sub-concentric
  common_scale_factors <- generate_scale_factors_all(k)
  common_scaled_df <- common_scale_factors %>%
    map_dfr(~create_scaled_df_sub(.x, df_circle, k))
  
  # For sub-concentric
  sub_concentric <- common_scaled_df %>%
    select(x_section = x_scaled, y_section = y_scaled, segment_id, sample_id, area_id = concentric_id) 
  
  # For super-sector
  super_sector <- create_sector_df(df_circle)
  
  # For super-concentric
  super_scale_factors <- generate_scale_factors_all_outside(k)
  super_scaled_df <- super_scale_factors %>%
    map_dfr(~create_scaled_df_super(.x, df_circle, k))
  
  super_concentric <- super_scaled_df %>%
    select(x_section = x_scaled, y_section = y_scaled, segment_id, sample_id, area_id = concentric_id) 
  
  # Modify the area_id column for each table
  sub_concentric <- sub_concentric %>% mutate(area_id = sapply(area_id, modify_area_id))
  super_concentric <- super_concentric %>% mutate(area_id = sapply(area_id, modify_area_id))
  
  return(list(sub_sector = sub_sector, sub_concentric = sub_concentric, super_sector = super_sector, super_concentric = super_concentric))
}