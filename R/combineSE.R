#' Combine SummarizedExperiment Objects
#'
#' @description
#' This function combines multiple SummarizedExperiment objects into a single
#' SummarizedExperiment object. It aligns features (rows) across all objects
#' based on the object with the most features, fills missing feature data with NA,
#' and concatenates colData from all objects.
#'
#' @param se_list A list of SummarizedExperiment objects to be combined.
#'        All objects should have an assay named "counts". The function
#'        assumes that all colData DataFrames can be combined with rbind,
#'        and if they have different columns, they need to be harmonized beforehand.
#'
#' @return A SummarizedExperiment object that represents the combination of all
#'         input objects. The returned object contains a unified feature set
#'         from the object with the most features, combined assay data, and
#'         concatenated colData reflecting metadata for each sample across
#'         the combined dataset.
#'
#' @export
#'
#' @examples
#' # Assuming you have a list of SummarizedExperiment objects: se_list
#' # combined_se <- combineSE(se_list)
#'
combineSE <- function(se_list) {
  # Step 1: Identify the SE object with the most rows (features)
  largest_idx <- which.max(sapply(se_list, function(x) nrow(rowData(x))))
  largest_se <- se_list[[largest_idx]]

  # Step 2: Prepare for combining by creating a matrix for the largest SE
  largest_mat <- assays(largest_se)[["counts"]]
  combined_mat <- largest_mat  # Start with the largest matrix

  # Initialize an empty list to store colData from each SE
  all_colData <- list()

  # Collect colData from the largest SE first
  all_colData[[largest_idx]] <- colData(largest_se)

  # Step 3: Loop through the other SE objects and align/combine them
  for (i in seq_along(se_list)) {
    if (i != largest_idx) {  # Skip the largest since it's already included
      current_se <- se_list[[i]]
      current_mat <- matrix(NA, nrow = nrow(largest_se), ncol = ncol(assays(current_se)[["counts"]]),
                            dimnames = list(rownames(largest_se), colnames(assays(current_se)[["counts"]])))

      # Find common features and fill the current matrix with data
      common_features <- intersect(rownames(largest_se), rownames(current_se))
      current_mat[common_features, ] <- assays(current_se)[["counts"]][common_features, ]

      # Combine with the existing matrix
      combined_mat <- cbind(combined_mat, current_mat)

      # Collect colData from the current SE
      all_colData[[i]] <- colData(current_se)
    }
  }

  # Combine all_colData into one DataFrame. This assumes each colData has the same columns.
  combined_colData <- do.call(rbind, all_colData)

  # Make sure the rownames of combined_colData match the colnames of the combined_mat
  rownames(combined_colData) <- colnames(combined_mat)

  # Create the combined SummarizedExperiment object with updated colData
  combined_se <- SummarizedExperiment(assays = list(counts = combined_mat), rowData = rowData(largest_se), colData = combined_colData)

  return(combined_se)
}
