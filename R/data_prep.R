# ==============================================================================
# File: R/data_prep.R
# Purpose: Clean raw field CSV data, standardise species/observer typos and
#          encoding glitches, sanitise GPS coordinates, validate angles,
#          compute perpendicular distances and transect effort, and attach
#          data-quality diagnostics.
# ==============================================================================

library(tidyverse)
library(sf)
library(Distance)

# ------------------------------------------------------------------------------
# Helper function: clean text labels and merge minor spelling mistakes
# ------------------------------------------------------------------------------
clean_labels <- function(x, default_label = "Unknown") {
  x_clean <- x |>
    as.character() |>
    str_remove_all("[[:cntrl:]]") |>         # Remove hidden control characters (e.g. \b)
    str_replace_all("C[\\(\\)]", "e") |>     # Fix UTF-8/Latin1 accent glitches (liC(ge, MottC)e)
    str_replace_all("[,;_-]+", " ") |>       # Standardise separators
    str_squish() |>
    str_to_sentence()

  x_clean <- if_else(is.na(x_clean) | x_clean == "", default_label, x_clean)

  freqs <- sort(table(x_clean), decreasing = TRUE)
  unique_vals <- names(freqs)

  if (length(unique_vals) <= 1) {
    return(x_clean)
  }

  canonical <- unique_vals
  for (i in 2:length(unique_vals)) {
    word <- unique_vals[i]
    max_edits <- case_when(
      nchar(word) >= 8 ~ 2L,
      nchar(word) >= 4 ~ 1L,
      TRUE ~ 0L
    )

    if (max_edits > 0L) {
      distances <- as.numeric(adist(word, canonical[1:(i - 1)], ignore.case = TRUE))
      best_match <- which.min(distances)
      if (distances[best_match] <= max_edits) {
        canonical[i] <- canonical[best_match]
      }
    }
  }

  canonical[match(x_clean, unique_vals)]
}

# ------------------------------------------------------------------------------
# Helper function: parse dirty numeric strings (e.g. "46..147" or "46,147")
# ------------------------------------------------------------------------------
parse_clean_numeric <- function(x) {
  x |>
    as.character() |>
    str_squish() |>
    str_replace_all(",", ".") |>
    str_replace_all("\\.{2,}", ".") |>
    as.numeric() |>
    suppressWarnings()
}

# ------------------------------------------------------------------------------
# Main function: prepare a Distance::ds() compatible data.frame + diagnostics
# ------------------------------------------------------------------------------
prepare_ds_data <- function(
  df,
  h_obs_m = 1.65,
  default_effort_km = 0.05,
  exclude_test_transects = TRUE,
  region_label = "Study_Area",
  area_km2 = 0.01
) {
  # Ensure optional columns exist even in minimal field CSV files
  if (!"species_common_name" %in% names(df)) df$species_common_name <- NA_character_
  if (!"species_code" %in% names(df)) df$species_code <- NA_character_
  if (!"observer" %in% names(df)) df$observer <- NA_character_
  if (!"count" %in% names(df)) df$count <- 1
  if (!"beaufort_id" %in% names(df)) df$beaufort_id <- NA_character_
  if (!"visibility" %in% names(df)) df$visibility <- NA_character_
  if (!"observation_timestamp" %in% names(df)) df$observation_timestamp <- seq_len(nrow(df))

  # Step 1: Standardise transect IDs, observer names, and GPS coordinates
  df <- df |>
    mutate(
      transect_name = str_squish(as.character(transect_name)),
      is_test_transect = str_detect(transect_name, regex("^test\\b", ignore_case = TRUE)),
      observer_clean = clean_labels(observer, default_label = "Unknown Observer") |>
        str_to_title(),
      lat_clean = parse_clean_numeric(gps_latitude),
      lon_clean = parse_clean_numeric(gps_longitude)
    )

  excluded_test_names <- df |>
    filter(!is.na(transect_name), is_test_transect) |>
    distinct(transect_name) |>
    pull(transect_name)

  if (isTRUE(exclude_test_transects)) {
    df <- df |>
      filter(is.na(transect_name) | !is_test_transect)
  }

  # Step 2: Filter USER detections and audit angle validity (degree)
  user_raw <- df |>
    filter(observation_type == "USER", !is.na(transect_name), transect_name != "") |>
    mutate(
      degree_num = parse_clean_numeric(degree),
      angle_status = case_when(
        is.na(degree_num) ~ "na",
        degree_num <= 0 ~ "le_0",
        degree_num > 90 ~ "gt_90",
        degree_num == 90 ~ "eq_90",
        TRUE ~ "valid"
      )
    )

  detections <- user_raw |>
    filter(angle_status %in% c("valid", "eq_90")) |>
    mutate(
      sp_name_clean = na_if(str_squish(as.character(species_common_name)), ""),
      sp_code_clean = na_if(str_squish(as.character(species_code)), ""),
      # Enrich generic "Bouchon" with species_code when common name lacks colour
      raw_species = case_when(
        !is.na(sp_name_clean) & str_to_lower(sp_name_clean) == "bouchon" & !is.na(sp_code_clean) ~
          paste("Bouchon", str_to_lower(sp_code_clean)),
        !is.na(sp_name_clean) ~ sp_name_clean,
        TRUE ~ sp_code_clean
      ),
      species = clean_labels(raw_species, default_label = "Unknown Species"),
      size = coalesce(parse_clean_numeric(count), 1),
      # Perpendicular distance in metres (degree = 90 gives d = 0 m on the line)
      distance = h_obs_m * tan((90 - degree_num) * pi / 180)
    ) |>
    filter(!is.na(distance), is.finite(distance), distance >= 0)

  # Step 3: Calculate transect effort (in kilometres) from sanitised GPS track points
  all_transects <- df |>
    distinct(transect_name) |>
    filter(!is.na(transect_name), transect_name != "")

  gps_raw <- df |>
    filter(observation_type == "GPS", !is.na(transect_name), transect_name != "")

  n_bad_gps <- sum(
    is.na(gps_raw$lat_clean) | is.na(gps_raw$lon_clean) |
    gps_raw$lat_clean < -90 | gps_raw$lat_clean > 90 |
    gps_raw$lon_clean < -180 | gps_raw$lon_clean > 180
  )

  gps_clean_pts <- gps_raw |>
    filter(
      !is.na(lat_clean), !is.na(lon_clean),
      between(lat_clean, -90, 90),
      between(lon_clean, -180, 180)
    ) |>
    arrange(transect_name, observation_timestamp) |>
    group_by(transect_name) |>
    filter(
      lat_clean != lag(lat_clean, default = first(lat_clean) + 999) |
      lon_clean != lag(lon_clean, default = first(lon_clean) + 999)
    ) |>
    filter(n() >= 2)

  if (nrow(gps_clean_pts) > 0) {
    effort_gps <- gps_clean_pts |>
      st_as_sf(coords = c("lon_clean", "lat_clean"), crs = 4326) |>
      summarise(do_union = FALSE, .groups = "drop") |>
      st_cast("LINESTRING") |>
      mutate(effort_km_gps = as.numeric(st_length(geometry)) / 1000) |>
      st_drop_geometry()
  } else {
    effort_gps <- tibble(
      transect_name = character(0),
      effort_km_gps = numeric(0)
    )
  }

  effort_by_transect <- all_transects |>
    left_join(effort_gps, by = "transect_name") |>
    mutate(
      use_fallback = is.na(effort_km_gps) | !is.finite(effort_km_gps) | effort_km_gps <= 0,
      Effort = if_else(use_fallback, default_effort_km, effort_km_gps),
      effort_source = if_else(use_fallback, "Fallback (missing/zero GPS)", "GPS")
    ) |>
    select(transect_name, Effort, effort_source)

  # Step 4: Assemble final data.frame in Distance::ds() format
  ds_data <- effort_by_transect |>
    left_join(
      detections |>
        select(
          transect_name,
          degree = degree_num,
          distance,
          size,
          species,
          obs_name = observer_clean,
          beaufort_id,
          visibility
        ),
      by = "transect_name"
    ) |>
    mutate(
      Region.Label = region_label,
      Area = area_km2,
      Sample.Label = transect_name,
      object = if_else(!is.na(distance), row_number(), NA_integer_)
    ) |>
    select(
      Region.Label,
      Area,
      Sample.Label,
      Effort,
      object,
      degree,
      distance,
      size,
      species,
      obs_name,
      beaufort_id,
      visibility,
      effort_source
    ) |>
    as.data.frame()

  # Attach quality audit metadata for Shiny diagnostic banners
  attr(ds_data, "diagnostics") <- list(
    n_user_raw = nrow(user_raw),
    n_valid = nrow(detections),
    n_na_degree = sum(user_raw$angle_status == "na"),
    n_le_0 = sum(user_raw$angle_status == "le_0"),
    n_gt_90 = sum(user_raw$angle_status == "gt_90"),
    n_eq_90 = sum(user_raw$angle_status == "eq_90"),
    n_bad_gps = n_bad_gps,
    excluded_test_names = if (isTRUE(exclude_test_transects)) excluded_test_names else character(0),
    n_transects = nrow(effort_by_transect),
    fallback_transects = effort_by_transect |>
      filter(effort_source != "GPS") |>
      pull(transect_name)
  )

  return(ds_data)
}