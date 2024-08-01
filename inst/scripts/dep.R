
#' Create sectors from a given data frame
#'
#' @param df A data frame to create sectors from.
#' @return A data frame of sectors.
create_sectors <- function(df) {
  df <- df[order(df$angle),]
  sectors <- list()

  for (i in 1:nrow(df)) {
    p1 <- df[i,]
    p2 <- if (i == nrow(df)) df[1,] else df[i+1,]

    sector <- data.frame(
      x_location_sector = c(p1$x_location, p2$x_location, p1$x_central, p1$x_location),
      y_location_sector = c(p1$y_location, p2$y_location, p1$y_central, p1$y_location),
      segment_id = p1$segment_id,
      sample_id = p1$sample_id,
      sector_id = paste0(p1$segment_id, "_", sprintf("%02d", i))
    )
    sectors[[i]] <- sector
  }

  do.call(rbind, sectors)
}


#' Create a scaled data frame
#'
#' @param s A numeric scaling factor.
#' @param df_circle A data frame to scale.
#' @param k A numeric value, defaulting to 8.
#' @return A scaled data frame.
create_scaled_df_sub <- function(s, df_circle, k = 5) {
  scaling_factor <- s * (k + 1)

  df_circle %>%
    group_by(segment_id) %>%
    mutate(
      x_diff = x_location - x_central,
      y_diff = y_location - y_central,
      x_scaled = x_central + s * x_diff,
      y_scaled = y_central + s * y_diff,
      scale = s,
      concentric_id = paste0(segment_id, "_", scaling_factor)
    ) %>%
    select(-x_diff, -y_diff) # Remove temporary columns
}


#' Create a scaled data frame
#'
#' @param s A numeric scaling factor.
#' @param df_circle A data frame to scale.
#' @param k A numeric value, defaulting to 8.
#' @return A scaled data frame.
create_scaled_df_super <- function(s, df_circle, k = 5) {
  scaling_factor <- s * (k + 1)-k

  df_circle %>%
    group_by(segment_id) %>%
    mutate(
      x_diff = x_location - x_central,
      y_diff = y_location - y_central,
      x_scaled = x_central + s * x_diff,
      y_scaled = y_central + s * y_diff,
      scale = s,
      concentric_id = paste0(segment_id, "_", scaling_factor)
    ) %>%
    select(-x_diff, -y_diff) # Remove temporary columns
}
