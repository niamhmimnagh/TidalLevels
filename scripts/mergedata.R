# ==============================================================================
# merge-mi-opw.R
#
# Creates a unified MI + OPW tide-gauge dataset.
#
# Tasks performed:
#   1. Reads the MI data, quality-controlled OPW data and OPW coordinates.
#   2. Creates a vertical-datum audit table.
#   3. Converts OPW data to the MI schema.
#   4. Sanitises station IDs.
#   5. Resolves the three shared gauges explicitly.
#   6. Writes one master merged Arrow dataset.
#   7. Writes the two agreed QC-analysis datasets.
#   8. Updates templates/stations.txt.
#
# Shared-gauge rule used:
#   - Ballycotton: use MI only.
#   - Dunmore East: use MI only.
#   - Malin: use OPW before the MI record begins, then MI.
#
# The original QC flags are retained in the master dataset.
# ==============================================================================

library(arrow)
library(dplyr)
library(stringr)
library(readr)
library(lubridate)
library(purrr)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

mi_file <- file.path(
  "data",
  "IrishNationalTideGaugeNetwork_712b_5843_9069_2025-10-17_st-id.parquet"
)

opw_qc_file <- file.path(
  "data",
  "OPW_Gauger_Data",
  "opw_gauger_qc.parquet"
)

opw_coords_file <- file.path(
  "data",
  "opw-station-coords.csv"
)

output_dir <- file.path(
  "data",
  "merged_tide_gauge_data"
)

master_dir <- file.path(
  output_dir,
  "mi_opw_master"
)

iteration_a_dir <- file.path(
  output_dir,
  "iteration_a_qc1"
)

iteration_b_dir <- file.path(
  output_dir,
  "iteration_b_qc01"
)

audit_dir <- file.path(
  output_dir,
  "audit"
)

stations_file <- file.path(
  "templates",
  "stations.txt"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Stop rather than silently overwrite an earlier merged dataset
stopifnot(
  "MI parquet file was not found." =
    file.exists(mi_file),
  
  "Quality-controlled OPW parquet file was not found." =
    file.exists(opw_qc_file),
  
  "OPW station-coordinate CSV was not found." =
    file.exists(opw_coords_file)
)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

# Produce station IDs that are safe for file names, Quarto files and RDS files.
#
# For example:
#   "St.Patrick's Quay (East)" -> "StPatricksQuayEast"
#   "Malin Head - Portmore Pier" -> "MalinHead-PortmorePier"
sanitize_st_id <- function(x) {
  
  x |>
    iconv(
      from = "",
      to = "ASCII//TRANSLIT"
    ) |>
    str_replace_all(
      "&",
      "And"
    ) |>
    str_replace_all(
      "[^A-Za-z0-9_-]",
      ""
    ) |>
    str_replace_all(
      "_+",
      "_"
    ) |>
    str_replace_all(
      "-+",
      "-"
    ) |>
    str_remove(
      "^[-_]+"
    ) |>
    str_remove(
      "[-_]+$"
    )
}


# Return the first column name found among several possible alternatives.
find_column <- function(data_names, candidates) {
  
  found <- candidates[
    candidates %in% data_names
  ]
  
  if (length(found) == 0L) {
    return(NA_character_)
  }
  
  found[[1]]
}


# Remove an existing Arrow dataset before rewriting it.
remove_existing_dataset <- function(path) {
  
  if (dir.exists(path)) {
    unlink(
      path,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  invisible(path)
}

# ------------------------------------------------------------------------------
# 3. Open source datasets
# ------------------------------------------------------------------------------

mi <- open_dataset(
  mi_file,
  format = "parquet"
)

opw_qc <- open_dataset(
  opw_qc_file,
  format = "parquet"
)

opw_coords_raw <- read_csv(
  opw_coords_file,
  show_col_types = FALSE
)

message(
  "MI rows:  ",
  format(mi$num_rows, big.mark = ",")
)

message(
  "OPW rows: ",
  format(opw_qc$num_rows, big.mark = ",")
)

# ------------------------------------------------------------------------------
# 4. Standardise the OPW coordinate table
# ------------------------------------------------------------------------------

coord_names <- names(opw_coords_raw)

coord_code_col <- find_column(
  coord_names,
  c(
    "code",
    "station_code",
    "station_ref",
    "station_reference",
    "hydrometric_station_reference"
  )
)

coord_station_col <- find_column(
  coord_names,
  c(
    "station",
    "station_name",
    "name"
  )
)

coord_lat_col <- find_column(
  coord_names,
  c(
    "latitude",
    "lat"
  )
)

coord_lon_col <- find_column(
  coord_names,
  c(
    "longitude",
    "lon",
    "long"
  )
)

if (
  is.na(coord_lat_col) ||
  is.na(coord_lon_col)
) {
  stop(
    paste(
      "Could not identify latitude and longitude columns in",
      opw_coords_file
    )
  )
}

if (
  is.na(coord_code_col) &&
  is.na(coord_station_col)
) {
  stop(
    paste(
      "The coordinate table must contain either an OPW station code",
      "or a station-name column."
    )
  )
}

opw_coords <- opw_coords_raw |>
  transmute(
    code = if (!is.na(coord_code_col)) {
      as.character(.data[[coord_code_col]])
    } else {
      NA_character_
    },
    
    station = if (!is.na(coord_station_col)) {
      as.character(.data[[coord_station_col]])
    } else {
      NA_character_
    },
    
    latitude = as.numeric(.data[[coord_lat_col]]),
    longitude = as.numeric(.data[[coord_lon_col]])
  ) |>
  filter(
    !is.na(latitude),
    !is.na(longitude)
  ) |>
  distinct()

# Make sure each OPW code has only one coordinate pair.
coordinate_check <- opw_coords |>
  filter(!is.na(code)) |>
  distinct(
    code,
    latitude,
    longitude
  ) |>
  count(
    code,
    name = "n_coordinate_pairs"
  )

bad_coordinates <- coordinate_check |>
  filter(n_coordinate_pairs != 1L)

if (nrow(bad_coordinates) > 0L) {
  
  print(bad_coordinates)
  
  stop(
    "At least one OPW code has more than one coordinate pair."
  )
}

write_csv(
  opw_coords,
  file.path(
    audit_dir,
    "opw_coordinates_used.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. Define the shared-gauge correspondence explicitly
# ------------------------------------------------------------------------------

shared_station_map <- tribble(
  ~opw_station,    ~mi_station,                    ~merge_rule,
  "Ballycotton",   "Ballycotton Harbour",          "MI only",
  "Dunmore East",  "Dunmore East Harbour",         "MI only",
  "Malin",         "Malin Head - Portmore Pier",   "OPW before MI; MI thereafter"
) |>
  mutate(
    canonical_station_id = mi_station,
    canonical_st_id = sanitize_st_id(mi_station)
  )

write_csv(
  shared_station_map,
  file.path(
    audit_dir,
    "shared_station_merge_rules.csv"
  )
)

# ------------------------------------------------------------------------------
# 6. Obtain station-level information needed for the datum audit
# ------------------------------------------------------------------------------

mi_station_summary <- mi |>
  group_by(
    station_id
  ) |>
  summarise(
    first_mi_time = min(
      time,
      na.rm = TRUE
    ),
    last_mi_time = max(
      time,
      na.rm = TRUE
    ),
    mi_latitude = mean(
      latitude,
      na.rm = TRUE
    ),
    mi_longitude = mean(
      longitude,
      na.rm = TRUE
    )
  ) |>
  collect()

opw_station_summary <- opw_qc |>
  group_by(
    code,
    station
  ) |>
  summarise(
    first_opw_time = min(
      datetime,
      na.rm = TRUE
    ),
    last_opw_time = max(
      datetime,
      na.rm = TRUE
    )
  ) |>
  collect()

# ------------------------------------------------------------------------------
# 7. Vertical-datum register
# ------------------------------------------------------------------------------

# The measurements support OD Malin consistency for the three shared gauges.
# They cannot independently establish the official datum of the new stations.
#
# Consequently:
#   - shared stations are marked as empirically consistent with OD Malin;
#   - new stations remain pending formal OPW confirmation.

datum_register <- opw_station_summary |>
  left_join(
    shared_station_map,
    by = c(
      "station" = "opw_station"
    )
  ) |>
  mutate(
    canonical_station_id = coalesce(
      canonical_station_id,
      station
    ),
    
    canonical_st_id = coalesce(
      canonical_st_id,
      sanitize_st_id(canonical_station_id)
    ),
    
    assumed_vertical_datum = "OD Malin",
    
    datum_status = case_when(
      !is.na(mi_station) ~
        "Empirically consistent with OD Malin at shared gauge",
      
      TRUE ~
        "Pending formal OPW confirmation"
    ),
    
    evidence = case_when(
      station == "Dunmore East" ~
        paste(
          "Shared MI/OPW series; approximately -0.04 m",
          "median/mean MI-minus-OPW offset"
        ),
      
      station == "Ballycotton" ~
        paste(
          "Shared MI/OPW series; approximately -0.04 m",
          "MI-minus-OPW offset after QC"
        ),
      
      station == "Malin" ~
        paste(
          "Shared MI/OPW series; approximately +0.01 m",
          "MI-minus-OPW offset"
        ),
      
      TRUE ~
        paste(
          "No MI series at this station;",
          "official datum confirmation required from OPW"
        )
    ),
    
    action_required = case_when(
      !is.na(mi_station) ~
        "Retain evidence and obtain written OPW confirmation if available",
      
      TRUE ~
        paste(
          "Confirm station datum and obtain a reduction constant",
          "if it is not OD Malin"
        )
    )
  ) |>
  select(
    code,
    opw_station = station,
    canonical_station_id,
    canonical_st_id,
    assumed_vertical_datum,
    datum_status,
    evidence,
    action_required,
    first_opw_time,
    last_opw_time
  ) |>
  arrange(
    datum_status,
    canonical_station_id
  )

write_csv(
  datum_register,
  file.path(
    audit_dir,
    "opw_vertical_datum_register.csv"
  )
)

print(
  datum_register |>
    count(
      datum_status
    )
)

# The merge can proceed provisionally, but produce a clear warning.
n_unconfirmed <- datum_register |>
  filter(
    datum_status == "Pending formal OPW confirmation"
  ) |>
  nrow()

if (n_unconfirmed > 0L) {
  warning(
    n_unconfirmed,
    paste(
      "OPW stations still require formal datum confirmation.",
      "See opw_vertical_datum_register.csv."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 8. Determine the start of the MI Malin record
# ------------------------------------------------------------------------------

malin_mi_start <- mi_station_summary |>
  filter(
    station_id == "Malin Head - Portmore Pier"
  ) |>
  pull(
    first_mi_time
  )

if (length(malin_mi_start) != 1L) {
  stop(
    "Could not identify exactly one MI start time for Malin."
  )
}

message(
  "The MI Malin series begins at: ",
  malin_mi_start
)

# ------------------------------------------------------------------------------
# 9. Join coordinates onto OPW records
# ------------------------------------------------------------------------------

opw_names <- names(opw_qc$schema)

opw_has_lat <- "latitude" %in% opw_names
opw_has_lon <- "longitude" %in% opw_names

# Use the coordinates already present in the QC parquet when available.
# Otherwise, join the coordinate table retrieved from waterlevel.ie.
if (opw_has_lat && opw_has_lon) {
  
  opw_with_coords <- opw_qc |>
    mutate(
      code = as.character(code)
    )
  
} else {
  
  # Prefer station code because it is the stable OPW identifier.
  if (!all(is.na(opw_coords$code))) {
    
    opw_coords_for_join <- opw_coords |>
      filter(
        !is.na(code)
      ) |>
      select(
        code,
        latitude,
        longitude
      ) |>
      distinct(
        code,
        .keep_all = TRUE
      )
    
    opw_with_coords <- opw_qc |>
      mutate(
        code = as.character(code)
      ) |>
      left_join(
        opw_coords_for_join,
        by = "code"
      )
    
  } else {
    
    opw_coords_for_join <- opw_coords |>
      filter(
        !is.na(station)
      ) |>
      select(
        station,
        latitude,
        longitude
      ) |>
      distinct(
        station,
        .keep_all = TRUE
      )
    
    opw_with_coords <- opw_qc |>
      left_join(
        opw_coords_for_join,
        by = "station"
      )
  }
}

# Verify coordinate coverage without collecting the full OPW dataset.
opw_missing_coordinates <- opw_with_coords |>
  filter(
    is.na(latitude) |
      is.na(longitude)
  ) |>
  distinct(
    code,
    station
  ) |>
  collect()

if (nrow(opw_missing_coordinates) > 0L) {
  
  print(opw_missing_coordinates)
  
  stop(
    paste(
      "Some OPW stations still have no coordinates.",
      "The merge has been stopped."
    )
  )
}

# ------------------------------------------------------------------------------
# 10. Standardise OPW to the MI schema
# ------------------------------------------------------------------------------

# The OPW QC file was described as having qc_flag in {1, 9}.
# Accept QC_Flag as an alternative spelling to make the code more robust.

if ("qc_flag" %in% opw_names) {
  
  opw_qc_column <- "qc_flag"
  
} else if ("QC_Flag" %in% opw_names) {
  
  opw_qc_column <- "QC_Flag"
  
} else {
  
  stop(
    paste(
      "Neither qc_flag nor QC_Flag was found in",
      opw_qc_file
    )
  )
}

opw_standardised <- opw_with_coords |>
  collect() |>
  mutate(
    station_id = case_when(
      station == "Ballycotton" ~
        "Ballycotton Harbour",
      
      station == "Dunmore East" ~
        "Dunmore East Harbour",
      
      station == "Malin" ~
        "Malin Head - Portmore Pier",
      
      TRUE ~
        station
    ),
    
    st_id = sanitize_st_id(
      station_id
    ),
    
    time = datetime,
    
    Water_Level_OD_Malin = value,
    
    QC_Flag = as.integer(
      .data[[opw_qc_column]]
    ),
    
    source = "OPW",
    
    source_station_id = station,
    
    source_station_code = as.character(code),
    
    datum = "OD Malin",
    
    datum_status = case_when(
      station %in% shared_station_map$opw_station ~
        "Empirically consistent with OD Malin",
      
      TRUE ~
        "Pending formal OPW confirmation"
    )
  ) |>
  select(
    time,
    station_id,
    st_id,
    latitude,
    longitude,
    Water_Level_OD_Malin,
    QC_Flag,
    source,
    source_station_id,
    source_station_code,
    datum,
    datum_status
  )
# ------------------------------------------------------------------------------
# 11. Standardise MI and repair st_id consistently
# ------------------------------------------------------------------------------

# Recalculate st_id for both sources with the same sanitisation function.
# This prevents punctuation from entering generated QMD and RDS file names.

mi_standardised <- mi |>
  collect() |>
  mutate(
    st_id = sanitize_st_id(station_id),
    source = "MI",
    source_station_id = station_id,
    source_station_code = NA_character_,
    datum = "OD Malin",
    datum_status = "Confirmed MI field"
  ) |>
  select(
    time,
    station_id,
    st_id,
    latitude,
    longitude,
    Water_Level_OD_Malin,
    QC_Flag,
    source,
    source_station_id,
    source_station_code,
    datum,
    datum_status
  )

# ------------------------------------------------------------------------------
# 11a. Enforce one canonical coordinate pair per merged station
# ------------------------------------------------------------------------------

# MI coordinates are preferred for stations represented by MI. OPW coordinates
# are used for OPW-only stations. Within a source, the most frequently occurring
# coordinate pair is selected.

mi_coordinate_candidates <- mi_standardised |>
  filter(
    !is.na(latitude),
    !is.na(longitude)
  ) |>
  group_by(
    station_id,
    st_id,
    latitude,
    longitude
  ) |>
  summarise(
    n_records = n()
  ) |>
  collect() |>
  mutate(
    coordinate_source = "MI",
    source_priority = 1L
  )

opw_coordinate_candidates <- opw_standardised |>
  filter(
    !is.na(latitude),
    !is.na(longitude)
  ) |>
  group_by(
    station_id,
    st_id,
    latitude,
    longitude
  ) |>
  summarise(
    n_records = n()
  ) |>
  collect() |>
  mutate(
    coordinate_source = "OPW",
    source_priority = 2L
  )

coordinate_candidates <- bind_rows(
  mi_coordinate_candidates,
  opw_coordinate_candidates
) |>
  arrange(
    station_id,
    source_priority,
    desc(n_records),
    latitude,
    longitude
  )

canonical_station_coordinates <- coordinate_candidates |>
  group_by(
    station_id,
    st_id
  ) |>
  slice_head(
    n = 1
  ) |>
  ungroup() |>
  select(
    station_id,
    st_id,
    latitude,
    longitude,
    coordinate_source,
    coordinate_n_records = n_records
  ) |>
  arrange(
    station_id
  )

## Manually verified corrections to source station coordinates.
## These are applied after selecting the canonical MI/OPW coordinate so that
## known source-coordinate errors do not propagate into the spatial analysis.
## Updated station coordinates supplied for the merged MI/OPW dataset
coordinate_corrections <- tribble(
  ~st_id,                    ~corrected_latitude, ~corrected_longitude,
  "AdelphiQuay",               52.2597,   -7.1024,
  "AranmoreIsland-Leabgarrow", 54.9894,   -8.4956,
  "BallycottonHarbour",        51.8281,   -8.0008,
  "BallyglassHarbour",         54.2533,   -9.8910,
  "Carlingford",               54.0523,   -6.1911,
  "CarrigaholtPier",           52.6005,   -9.6982,
  "CastletownberePort",        51.6493,   -9.9031,
  "CurrachClub",               51.9015,   -8.4360,
  "DingleHarbour",             52.1385,  -10.2774,
  "DunmoreEastHarbour",        52.1475,   -6.9908,
  "EnniscronePier",            54.2202,   -9.0985,
  "Fenit",                     52.2697,   -9.8521,
  "Foynes",                    52.6130,   -9.1055,
  "GalwayPort",                53.2691,   -9.0482,
  "HowthWaterLevel1",          53.3918,   -6.0683,
  "HowthWaterLevel2",          53.3918,   -6.0683,
  "Inishmore",                 53.1178,   -9.6667,
  "Kenmare",                   51.8718,   -9.5884,
  "KillybegsPort",             54.6364,   -8.4396,
  "KilrushLough",              52.6316,   -9.5018,
  "LimerickDock",              52.6582,   -8.6438,
  "MalinHead-PortmorePier",    55.3716,   -7.3347,
  "MercyHospital",             51.90006,  -8.4838,
  "Pope'sQuay",                 51.9010,   -8.4735,
  "PortOriel",                 53.7990,   -6.2217,
  "RingaskiddyNMCI",           51.8349,   -8.3045,
  "RossaveelPier",             53.2669,   -9.5621,
  "Rosslare",                  52.2546,   -6.3349,
  "ShannonAirport",            52.6784,   -8.9173,
  "Sligo",                     54.3046,   -8.5688,
  "StPatrick'sQuay(East)",        51.9002,   -8.4639,
  "UnionHallHarbor",           51.5585,   -9.1329,
  "UnionHallHarbor2",          51.5585,   -9.1329,
  "WexfordHarbour",            52.3340,   -6.4533,
  "YoughalQuay",               51.9568,   -7.8468
)

## Check that every supplied station ID exists in the merged dataset
unmatched_corrections <- coordinate_corrections |>
  anti_join(canonical_station_coordinates, by = "st_id")

if (nrow(unmatched_corrections) > 0) {
  print(unmatched_corrections)
  warning("Some coordinate corrections could not be matched to station IDs.")
}

## Identify stations without a supplied coordinate correction
uncorrected_stations <- canonical_station_coordinates |>
  anti_join(coordinate_corrections, by = "st_id") |>
  select(station_id, st_id, latitude, longitude)

print(uncorrected_stations)

canonical_station_coordinates <- canonical_station_coordinates |>
  left_join(coordinate_corrections, by = "st_id") |>
  mutate(
    coordinate_corrected = !is.na(corrected_latitude),
    latitude_original = latitude,
    longitude_original = longitude,
    latitude = coalesce(corrected_latitude, latitude),
    longitude = coalesce(corrected_longitude, longitude),
    coordinate_source = if_else(coordinate_corrected, "Manually verified correction", coordinate_source)
  ) |>
  select(-corrected_latitude, -corrected_longitude)

if (nrow(canonical_station_coordinates) == 0L) {
  stop(
    "No valid canonical station coordinates were available."
  )
}

write_csv(
  coordinate_candidates,
  file.path(
    audit_dir,
    "station_coordinate_candidates.csv"
  )
)

write_csv(
  canonical_station_coordinates,
  file.path(
    audit_dir,
    "canonical_station_coordinates.csv"
  )
)

canonical_coordinate_join <- canonical_station_coordinates |>
  select(
    station_id,
    st_id,
    latitude,
    longitude
  )

mi_standardised <- mi_standardised |>
  select(
    -latitude,
    -longitude
  ) |>
  left_join(
    canonical_coordinate_join,
    by = c(
      "station_id",
      "st_id"
    )
  ) |>
  select(
    time,
    station_id,
    st_id,
    latitude,
    longitude,
    Water_Level_OD_Malin,
    QC_Flag,
    source,
    source_station_id,
    source_station_code,
    datum,
    datum_status
  )

opw_standardised <- opw_standardised |>
  select(
    -latitude,
    -longitude
  ) |>
  left_join(
    canonical_coordinate_join,
    by = c(
      "station_id",
      "st_id"
    )
  ) |>
  select(
    time,
    station_id,
    st_id,
    latitude,
    longitude,
    Water_Level_OD_Malin,
    QC_Flag,
    source,
    source_station_id,
    source_station_code,
    datum,
    datum_status
  )

# ------------------------------------------------------------------------------
# 12. Deduplicate the shared gauges
# ------------------------------------------------------------------------------

# Policy:
#
# Ballycotton:
#   Drop OPW entirely because MI is already quality controlled.
#
# Dunmore East:
#   Drop OPW entirely because MI is already quality controlled.
#
# Malin:
#   Keep OPW only before the first MI observation.
#   Keep MI from its normal starting date onwards.
#
# This prevents overlapping records from being double counted.

opw_deduplicated <- opw_standardised |>
  filter(
    !station_id %in% c(
      "Ballycotton Harbour",
      "Dunmore East Harbour"
    )
  ) |>
  filter(
    station_id != "Malin Head - Portmore Pier" |
      time < malin_mi_start
  )

deduplication_audit <- tibble(
  canonical_station_id = c(
    "Ballycotton Harbour",
    "Dunmore East Harbour",
    "Malin Head - Portmore Pier"
  ),
  
  mi_action = c(
    "Keep complete MI series",
    "Keep complete MI series",
    "Keep MI series from its original start onwards"
  ),
  
  opw_action = c(
    "Drop complete OPW series",
    "Drop complete OPW series",
    paste(
      "Keep OPW records before",
      format(malin_mi_start)
    )
  ),
  
  reason = c(
    "MI is quality controlled; avoid duplicate observations",
    "MI is quality controlled; avoid duplicate observations",
    "OPW extends the record earlier without overlapping the retained MI period"
  )
)

write_csv(
  deduplication_audit,
  file.path(
    audit_dir,
    "shared_gauge_deduplication_audit.csv"
  )
)

# ------------------------------------------------------------------------------
# 13. Bind MI and OPW
# ------------------------------------------------------------------------------

merged_master <- union_all(
  mi_standardised,
  opw_deduplicated
)


# Verify that each canonical station has exactly one coordinate pair.
merged_coordinate_check <- merged_master |>
  distinct(
    station_id,
    st_id,
    latitude,
    longitude
  ) |>
  collect() |>
  count(
    station_id,
    st_id,
    name = "n_coordinate_pairs"
  )

bad_merged_coordinates <- merged_coordinate_check |>
  filter(
    n_coordinate_pairs != 1L
  )

if (nrow(bad_merged_coordinates) > 0L) {
  print(
    bad_merged_coordinates
  )
  
  write_csv(
    bad_merged_coordinates,
    file.path(
      audit_dir,
      "bad_merged_station_coordinates.csv"
    )
  )
  
  stop(
    "At least one merged station still has multiple coordinate pairs."
  )
}

write_csv(
  canonical_station_coordinates,
  file.path(
    audit_dir,
    "merged_station_coordinates_final.csv"
  )
)

# Check for duplicated source-time records at a station.
#
# A duplicate should not normally exist after the shared-gauge rules.
# This only collects the resulting duplicate summary, not the full dataset.

duplicate_station_times <- merged_master |>
  count(
    station_id,
    time,
    name = "n_records"
  ) |>
  filter(
    n_records > 1L
  ) |>
  collect()

if (nrow(duplicate_station_times) > 0L) {
  
  write_csv(
    duplicate_station_times,
    file.path(
      audit_dir,
      "duplicate_station_times.csv"
    )
  )
  
  warning(
    nrow(duplicate_station_times),
    paste(
      "duplicated station-time combinations remain.",
      "See duplicate_station_times.csv."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 14. Create the two agreed QC iterations
# ------------------------------------------------------------------------------

# Iteration A:
#   MI QC_Flag = 1
#   OPW screened-good records, QC_Flag = 1
iteration_a <- merged_master |>
  filter(QC_Flag == 1L) |>
  filter(
    !(st_id == "Inishmore" & time < as.POSIXct("2016-01-01", tz = "UTC")),
    !(st_id == "Sligo" & time < as.POSIXct("2018-01-01", tz = "UTC")),
    !(st_id == "SkerriesHarbour" & time < as.POSIXct("2013-01-01", tz = "UTC"))
  )


# Iteration B:
#   MI QC_Flag in {0, 1}
#   OPW screened-good records only
#
# OPW records flagged 9 are excluded. The OPW QC process does not produce
# unchecked flag 0 records.
iteration_b <- merged_master |>
  filter(
    (
      source == "MI" &
        QC_Flag %in% c(
          0L,
          1L
        )
    ) |
      (
        source == "OPW" &
          QC_Flag == 1L
      )
  )

# ------------------------------------------------------------------------------
# 15. Write the Arrow datasets
# ------------------------------------------------------------------------------

# Arrow datasets are used instead of collecting approximately 39 million rows
# into memory. Each output is a directory containing parquet parts and can be
# opened using arrow::open_dataset().

walk(
  c(
    master_dir,
    iteration_a_dir,
    iteration_b_dir
  ),
  remove_existing_dataset
)

write_dataset(
  merged_master,
  path = master_dir,
  format = "parquet",
  partitioning = c(
    "source"
  ),
  existing_data_behavior = "overwrite"
)

write_dataset(
  iteration_a,
  path = iteration_a_dir,
  format = "parquet",
  partitioning = c(
    "source"
  ),
  existing_data_behavior = "overwrite"
)

write_dataset(
  iteration_b,
  path = iteration_b_dir,
  format = "parquet",
  partitioning = c(
    "source"
  ),
  existing_data_behavior = "overwrite"
)

# ------------------------------------------------------------------------------
# 16. Produce station and row-count audits
# ------------------------------------------------------------------------------

master_station_audit <- merged_master |>
  group_by(
    source,
    station_id,
    st_id
  ) |>
  summarise(
    first = min(
      time,
      na.rm = TRUE
    ),
    last = max(
      time,
      na.rm = TRUE
    ),
    n = n(),
    latitude = mean(
      latitude,
      na.rm = TRUE
    ),
    longitude = mean(
      longitude,
      na.rm = TRUE
    )
  ) |>
  collect() |>
  arrange(
    station_id,
    source,
    first
  )

write_csv(
  master_station_audit,
  file.path(
    audit_dir,
    "merged_station_audit.csv"
  )
)

iteration_counts <- bind_rows(
  iteration_a |>
    summarise(
      n = n()
    ) |>
    collect() |>
    mutate(
      iteration = "A: MI QC 1 plus OPW QC 1"
    ),
  
  iteration_b |>
    summarise(
      n = n()
    ) |>
    collect() |>
    mutate(
      iteration = "B: MI QC 0/1 plus OPW QC 1"
    )
) |>
  select(
    iteration,
    n
  )

write_csv(
  iteration_counts,
  file.path(
    audit_dir,
    "iteration_row_counts.csv"
  )
)

print(iteration_counts)

# ------------------------------------------------------------------------------
# 17. Check that every station ID is safe and unique
# ------------------------------------------------------------------------------

station_ids <- merged_master |>
  distinct(
    station_id,
    st_id
  ) |>
  collect() |>
  arrange(
    st_id
  )

unsafe_st_ids <- station_ids |>
  filter(
    is.na(st_id) |
      st_id == "" |
      str_detect(
        st_id,
        "[^A-Za-z0-9_-]"
      )
  )

if (nrow(unsafe_st_ids) > 0L) {
  
  print(unsafe_st_ids)
  
  stop(
    "At least one st_id remains unsafe for use in file names."
  )
}

duplicated_st_ids <- station_ids |>
  count(
    st_id,
    name = "n_station_names"
  ) |>
  filter(
    n_station_names > 1L
  )

if (nrow(duplicated_st_ids) > 0L) {
  
  print(
    station_ids |>
      semi_join(
        duplicated_st_ids,
        by = "st_id"
      )
  )
  
  stop(
    paste(
      "Different station names produced the same sanitised st_id.",
      "These IDs must be resolved manually."
    )
  )
}

write_csv(
  station_ids,
  file.path(
    audit_dir,
    "station_id_lookup.csv"
  )
)

# ------------------------------------------------------------------------------
# 18. Update templates/stations.txt
# ------------------------------------------------------------------------------

# Read existing stations when the file already exists.
existing_stations <- if (file.exists(stations_file)) {
  
  read_lines(
    stations_file
  ) |>
    str_trim() |>
    discard(
      ~ .x == ""
    )
  
} else {
  
  character()
}

all_station_ids <- sort(
  unique(
    c(
      existing_stations,
      station_ids$st_id
    )
  )
)

write_lines(
  all_station_ids,
  stations_file
)

message(
  "Updated ",
  stations_file,
  " with ",
  length(all_station_ids),
  " unique station IDs."
)

# ------------------------------------------------------------------------------
# 19. Final checks
# ------------------------------------------------------------------------------

master_check <- open_dataset(
  master_dir,
  format = "parquet"
)

iteration_a_check <- open_dataset(
  iteration_a_dir,
  format = "parquet"
)

iteration_b_check <- open_dataset(
  iteration_b_dir,
  format = "parquet"
)

final_check <- bind_rows(
  master_check |>
    summarise(
      n = n(),
      n_stations = n_distinct(st_id),
      first = min(time),
      last = max(time)
    ) |>
    collect() |>
    mutate(
      dataset = "Master"
    ),
  
  iteration_a_check |>
    summarise(
      n = n(),
      n_stations = n_distinct(st_id),
      first = min(time),
      last = max(time)
    ) |>
    collect() |>
    mutate(
      dataset = "Iteration A"
    ),
  
  iteration_b_check |>
    summarise(
      n = n(),
      n_stations = n_distinct(st_id),
      first = min(time),
      last = max(time)
    ) |>
    collect() |>
    mutate(
      dataset = "Iteration B"
    )
) |>
  select(
    dataset,
    everything()
  )

print(final_check)

write_csv(
  final_check,
  file.path(
    audit_dir,
    "final_dataset_check.csv"
  )
)

message(
  paste(
    "\nMerge completed successfully.",
    "\nMaster dataset:     ", master_dir,
    "\nIteration A:        ", iteration_a_dir,
    "\nIteration B:        ", iteration_b_dir,
    "\nAudit outputs:      ", audit_dir,
    "\nStation list:       ", stations_file
  )
)