## =============================================================================
## EXPORT MI/OPW TIDAL ANALYSIS RESULTS TO GEOPACKAGE
##
## Standalone script.
## Run from a clean R session after the final m1st prediction models have
## been generated using the corrected Ireland barrier and Shannon Estuary.
##
## Exports:
##
##   1. stations
##      One point per tide-gauge station.
##
##   2. tidal_features
##      Observed/derived station-level tidal characteristics.
##
##   3. pred_HAT
##   4. pred_MHWS
##   5. pred_MHWN
##   6. pred_MSL
##   7. pred_MLWN
##   8. pred_MLWS
##   9. pred_LAT
##
##      Each prediction layer contains:
##      year_idx, year_mid, longitude, latitude,
##      mean, sd, q025, q975 and geometry.
##
## Spatial setup:
##   - corrected Ireland barrier:
##       data/spatial/ireland_new_barrier.gpkg
##
##   - explicit Shannon Estuary water polygon:
##       data/spatial/shannon_estuary.gpkg
##
##   - general 2-km prediction grid over the marine domain
##
##
## Output:
##   outputs/gis/MI_OPW_tidal_analysis_v2.gpkg
## =============================================================================


## =============================================================================
## 1. PACKAGES
## =============================================================================

library(tidyverse)
library(sf)
library(INLA)
library(fmesher)
library(INLAspacetime)
library(stExtremes)
library(here)
itmproj <- st_crs(2157)
source(here::here("scripts", "_functions.R"))


## =============================================================================
## 2. SETTINGS
## =============================================================================

const_levels <- c("HAT", "MHWS", "MHWN", "MSL", "MLWN", "MLWS", "LAT")

prediction_model <- "m1st"

gis_dir <- here::here("outputs", "gis")
dir.create(gis_dir, recursive = TRUE, showWarnings = FALSE)

gpkg_file <- here::here(
  "outputs",
  "gis",
  "MI_OPW_tidal_analysis_v2.gpkg"
)

if (file.exists(gpkg_file)) {
  file.remove(gpkg_file)
}


## =============================================================================
## 3. READ SAVED SPATIO-TEMPORAL INPUT DATA
## =============================================================================

st_files <- list.files(
  path = here::here("rds"),
  pattern = "^inla-st-data-.*\\.rds$",
  full.names = TRUE
)

if (length(st_files) == 0) {
  stop(
    paste0(
      "No files matching 'inla-st-data-*.rds' were found in:\n",
      here::here("rds")
    )
  )
}

cat("\nReading", length(st_files), "spatio-temporal input files...\n")

da <- map_dfr(st_files, readRDS)

required_columns <- c(
  "st_id",
  "constituent",
  "year_mid",
  "value",
  "lat",
  "lon"
)

missing_columns <- setdiff(required_columns, names(da))

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}


## =============================================================================
## 4. CHECK TIDAL FEATURES
## =============================================================================

available_constituents <- unique(as.character(da$constituent))
missing_constituents <- setdiff(const_levels, available_constituents)

if (length(missing_constituents) > 0) {
  warning(
    paste0(
      "The following tidal features are missing: ",
      paste(missing_constituents, collapse = ", ")
    )
  )
}

da <- da |>
  filter(as.character(constituent) %in% const_levels) |>
  mutate(
    constituent = factor(
      as.character(constituent),
      levels = const_levels
    )
  )


## =============================================================================
## 5. RECONSTRUCT CORRECT STATION COORDINATES
## =============================================================================

station_coordinates <- da |>
  filter(!is.na(st_id), !is.na(lat), !is.na(lon)) |>
  distinct(st_id, lat, lon)

coordinate_check <- station_coordinates |>
  count(st_id, name = "n_coordinate_pairs")

if (any(coordinate_check$n_coordinate_pairs != 1L)) {
  print(
    coordinate_check |>
      filter(n_coordinate_pairs != 1L)
  )
  
  stop(
    "The source data contain multiple coordinate pairs for a station."
  )
}

cat(
  "Found",
  nrow(station_coordinates),
  "unique station locations.\n"
)


## =============================================================================
## 6. WRITE STATION LOCATIONS
## =============================================================================

stations_gis <- station_coordinates |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(itmproj)

station_xy <- st_coordinates(stations_gis)

stations_gis <- stations_gis |>
  mutate(Easting = station_xy[, "X"], Northing = station_xy[, "Y"])

st_write(
  stations_gis,
  gpkg_file,
  layer = "stations",
  quiet = TRUE
)

cat("Written layer: stations\n")


## =============================================================================
## 7. CREATE STATION-LEVEL TIDAL FEATURE LAYER
## =============================================================================

tidal_features_wide <- da |>
  filter(!is.na(value), !is.na(lat), !is.na(lon)) |>
  mutate(constituent = as.character(constituent)) |>
  select(st_id, year_mid, constituent, value, lat, lon) |>
  group_by(st_id, year_mid, constituent, lat, lon) |>
  summarise(value = first(value), .groups = "drop") |>
  pivot_wider(
    names_from = constituent,
    values_from = value
  )

stopifnot(!anyDuplicated(names(tidal_features_wide)))

tidal_features_gis <- tidal_features_wide |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(itmproj)

tidal_xy <- st_coordinates(tidal_features_gis)

tidal_features_gis <- tidal_features_gis |>
  mutate(Easting = tidal_xy[, "X"], Northing = tidal_xy[, "Y"]) |>
  rename(StationID = st_id, Year = year_mid)

st_write(
  tidal_features_gis,
  gpkg_file,
  layer = "tidal_features",
  append = FALSE,
  quiet = TRUE
)

cat("Written layer: tidal_features\n")


## =============================================================================
## 8. RECONSTRUCT DATA USED BY SPATIO-TEMPORAL MODELS
## =============================================================================

yrs <- range(da$year_mid, na.rm = TRUE)
yrs_seq <- seq(yrs[1], yrs[2])

db <- da |>
  split(~ constituent)

db <- db[
  intersect(const_levels, names(db))
]

db <- imap(
  db,
  \(x, nm) {
    x |>
      select(-lat, -lon) |>
      complete(
        year_mid = yrs_seq,
        st_id = st_id
      ) |>
      left_join(
        station_coordinates,
        by = "st_id"
      ) |>
      mutate(
        constituent = factor(
          nm,
          levels = const_levels
        )
      ) |>
      arrange(year_mid, st_id)
  }
)


## =============================================================================
## 9. CREATE TEMPORAL INDEX
## =============================================================================

db <- map(
  db,
  \(x) {
    x |>
      mutate(
        year_idx = as.numeric(as.factor(year_mid))
      )
  }
)


## =============================================================================
## 10. COORDINATE REFERENCE SYSTEMS
## =============================================================================

proj <- st_crs(
  "+proj=longlat +datum=WGS84"
)

kmproj <- st_crs(
  "+proj=utm +zone=30 +datum=WGS84 +units=km"
)



## =============================================================================
## 11. LOAD CORRECTED SPATIAL DOMAIN
## =============================================================================

data(area, package = "stExtremes")

barrier <- st_read(here::here("data", "spatial", "ireland_new_barrier.gpkg"), quiet = TRUE)
shannon <- st_read(here::here("data", "spatial", "shannon_estuary.gpkg"), quiet = TRUE)

area <- st_make_valid(st_as_sf(area))
barrier <- st_make_valid(barrier) |> st_transform(st_crs(area))
shannon <- st_make_valid(shannon) |> st_transform(st_crs(area))

## Extend the original marine domain into the Shannon Estuary
area <- st_union(area, shannon) |> st_make_valid()

area_km <- st_transform(area, kmproj)
barrier_km <- st_transform(barrier, kmproj)

## Exact marine model domain used in model fitting
area_model <- st_difference(
  st_make_valid(area_km),
  st_union(st_make_valid(barrier_km))
) |>
  st_make_valid()

area_buffered <- area_model


## =============================================================================
## 12. CREATE SPATIAL DATA OBJECTS
## =============================================================================

db_st <- map(
  db,
  \(x) {
    st_as_sf(
      x,
      coords = c("lon", "lat"),
      crs = proj
    )
  }
)

db_st_km <- map(
  db_st,
  \(x) {
    st_transform(x, kmproj)
  }
)

db_st <- map(
  db_st,
  add_coords_to_sf
)

db_st_km <- map(
  db_st_km,
  add_coords_to_sf
)


## =============================================================================
## 13. CENTER MODEL VARIABLES
## =============================================================================

scales_db_st <- map(
  db_st,
  \(x) {
    center_scale_vars(
      x,
      return_params = TRUE
    )
  }
)

db_st <- map2(
  db_st,
  scales_db_st,
  \(x, y) {
    bind_cols(
      x,
      y$data
    )
  }
)

scales_db_st_km <- map(
  db_st_km,
  \(x) {
    center_scale_vars(
      x,
      return_params = TRUE
    )
  }
)

db_st_km <- map2(
  db_st_km,
  scales_db_st_km,
  \(x, y) {
    bind_cols(
      x,
      y$data
    )
  }
)


## =============================================================================
## 14. RECONSTRUCT FINAL SPATIAL MESH
##
## IMPORTANT: this is copied from the modelling code and must remain identical
## to the mesh used when the saved prediction models were generated.
## =============================================================================

r2 <- matrix(st_bbox(area_km), nrow = 2) |> apply(1, diff)
rr <- mean(r2)

max.edge.original <- (rr / 15) * c(1, 5)
bound.outer <- max.edge.original[2] * 1.1
max.edge <- c(15, 70)

mesh0 <- fm_mesh_2d_inla(
  boundary = area_buffered,
  max.edge = max.edge[1] / 2,
  cutoff = max.edge[1] / 4
)

mesh0_xy <- mesh0$loc[, 1:2, drop = FALSE]

mesh <- fm_mesh_2d_inla(
  loc = mesh0_xy,
  max.edge = max.edge[2],
  offset = bound.outer,
  cutoff = max.edge[1] / 2
)

fm_crs(mesh) <- kmproj

cat("\nMesh vertices:", mesh$n, "\n")
cat("Mesh triangles:", nrow(mesh$graph$tv), "\n")


## =============================================================================
## 15. RECONSTRUCT BARRIER MODEL
## =============================================================================

tri_barrier <- tri_on_mesh(
  mesh,
  area_model,
  inverse = TRUE
)

prior_sigma <- c(1, 0.01)
prior_range <- c(50, 0.01)

barrier_model <- INLAspacetime::barrierModel.define(
  mesh,
  barrier.triangles = tri_barrier,
  prior.range = prior_range,
  prior.sigma = prior_sigma,
  constr = TRUE
)


## =============================================================================
## 16. RECONSTRUCT SPATIO-TEMPORAL PROJECTOR MATRICES
## =============================================================================

Ast.data <- map(
  db_st_km,
  \(x) {
    inla.spde.make.A(
      mesh,
      loc = st_coordinates(x),
      group = x$year_idx,
      n.group = length(unique(x$year_idx))
    )
  }
)

stopifnot(
  all(
    map2_lgl(
      Ast.data,
      db_st_km,
      \(A, x) {
        isTRUE(
          all.equal(
            nrow(x),
            sum(A)
          )
        )
      }
    )
  )
)


## =============================================================================
## 17. RECONSTRUCT LATENT SPATIO-TEMPORAL INDICES
## =============================================================================

w.index <- map(
  db_st_km,
  \(x) {
    create_indices(
      x$year_idx,
      mesh
    )
  }
)


## =============================================================================
## 18. RECONSTRUCT ESTIMATION STACKS
## =============================================================================

stk.st <- pmap(
  list(
    db_st_km,
    w.index,
    Ast.data
  ),
  \(x, w, A) {
    inla_stack_st(
      y = x$value,
      w = w,
      trend = x$trendc,
      A = A,
      tag = "est.st"
    )
  }
)

names(stk.st) <- names(db_st_km)


## =============================================================================
## 19. RECONSTRUCT PREDICTION GRID
##
## This is deliberately identical to the grid used in the modelling code.
## No separate local river or estuary grids are added here.
## =============================================================================

grid_pts <- st_make_grid(area_model, cellsize = 2, what = "centers")

id.grid.in <- unlist(st_contains(area_model, grid_pts))
n.pred <- length(id.grid.in)

grid_pts_inside <- grid_pts[id.grid.in]
pred.coords <- st_coordinates(grid_pts_inside)

pred_coord_df <- as.data.frame(pred.coords)
names(pred_coord_df) <- c("X", "Y")

cat("\nPrediction locations per year:", n.pred, "\n")


## =============================================================================
## 21. RECONSTRUCT TEMPORAL PREDICTION GRID
## =============================================================================

prediction_year_idx <- sort(unique(db_st_km[[1]]$year_idx))
n_years <- length(prediction_year_idx)

## Use the same function and ordering as the modelling/prediction code
expd_grid <- tidyr::expand_grid(pred_coord_df, year_idx = prediction_year_idx)

stopifnot(nrow(expd_grid) == n.pred * n_years)

cat("Prediction years:", n_years, "\n")
cat("Total space-time prediction rows:", nrow(expd_grid), "\n")


## =============================================================================
## 23. YEAR LOOKUP
## =============================================================================

year_lookup <- db_st_km[[1]] |>
  st_drop_geometry() |>
  distinct(year_idx, year_mid) |>
  arrange(year_idx)

stopifnot(
  nrow(year_lookup) ==
    n_years
)


## =============================================================================
## 24. RECONSTRUCT PREDICTION PROJECTOR MATRIX
## =============================================================================

pred_loc <- as.matrix(expd_grid[, c("X", "Y")])

A.pred <- inla.spde.make.A(
  mesh,
  loc = pred_loc,
  group = expd_grid$year_idx,
  n.group = n_years
)

stopifnot(isTRUE(all.equal(
  as.numeric(sum(A.pred)),
  as.numeric(n.pred * n_years)
)))

cat("Prediction projector matrix constructed successfully.\n")


## =============================================================================
## 25. READ SAVED PREDICTION MODELS AND EXTRACT RESULTS
##
## IMPORTANT:
## The prediction RDS files must have been generated using the same mesh,
## prediction grid, year ordering and stack construction above.
## =============================================================================

df_pred <- vector(
  mode = "list",
  length = length(names(db_st_km))
)

names(df_pred) <- names(
  db_st_km
)

for (const in names(db_st_km)) {
  
  cat(
    "\nExtracting",
    const,
    "predictions...\n"
  )
  
  prediction_file <- here::here(
    "rds",
    sprintf(
      "tide-models-st-PRED-FIXED-%s-%s.rds",
      const,
      prediction_model
    )
  )
  
  if (!file.exists(prediction_file)) {
    warning(
      sprintf(
        "Prediction model not found for %s: %s",
        const,
        prediction_file
      )
    )
    
    df_pred[[const]] <- NULL
    
    next
  }
  
  
  ## ---------------------------------------------------------------------------
  ## Prediction stack
  ## ---------------------------------------------------------------------------
  
  stk.pred.i <- inla.stack(
    data = list(
      y = NA
    ),
    effects = list(
      w.index[[const]],
      data.frame(
        m = rep(
          1,
          nrow(A.pred)
        ),
        trendc = expd_grid$year_idx -
          mean(expd_grid$year_idx)
      )
    ),
    A = list(
      A.pred,
      1
    ),
    tag = "pred"
  )
  
  stack.full.i <- inla.stack(
    stk.st[[const]],
    stk.pred.i
  )
  
  idx.pred.i <- inla.stack.index(
    stack.full.i,
    tag = "pred"
  )$data
  
  stopifnot(
    length(idx.pred.i) ==
      nrow(expd_grid)
  )
  
  
  ## ---------------------------------------------------------------------------
  ## Read saved prediction model
  ## ---------------------------------------------------------------------------
  
  pred_model_i <- readRDS(
    prediction_file
  )
  
  fitted_summary_i <- pred_model_i$summary.fitted.values
  
  if (max(idx.pred.i) > nrow(fitted_summary_i)) {
    stop(
      paste0(
        "\nThe reconstructed prediction stack does not match the saved ",
        "prediction model for ",
        const,
        ".\n\n",
        "Maximum reconstructed prediction index: ",
        max(idx.pred.i),
        "\n",
        "Rows in saved fitted-value summary: ",
        nrow(fitted_summary_i),
        "\n\n",
        "Do not export the GIS data until the saved prediction model and ",
        "the reconstructed prediction grid use exactly the same setup."
      )
    )
  }
  
  
  ## ---------------------------------------------------------------------------
  ## Extract posterior summaries
  ## ---------------------------------------------------------------------------
  
  pred_mean_i <- fitted_summary_i[
    idx.pred.i,
    "mean",
    drop = TRUE
  ]
  
  pred_sd_i <- fitted_summary_i[
    idx.pred.i,
    "sd",
    drop = TRUE
  ]
  
  pred_q025_i <- fitted_summary_i[
    idx.pred.i,
    "0.025quant",
    drop = TRUE
  ]
  
  pred_q975_i <- fitted_summary_i[
    idx.pred.i,
    "0.975quant",
    drop = TRUE
  ]
  
  stopifnot(
    length(pred_mean_i) == nrow(expd_grid),
    length(pred_sd_i) == nrow(expd_grid),
    length(pred_q025_i) == nrow(expd_grid),
    length(pred_q975_i) == nrow(expd_grid)
  )
  
  
  ## ---------------------------------------------------------------------------
  ## Attach predictions to their actual coordinates
  ## ---------------------------------------------------------------------------
  
  df_pred[[const]] <- tibble(
    X = expd_grid$X,
    Y = expd_grid$Y,
    year_idx = expd_grid$year_idx,
    mean = pred_mean_i,
    sd = pred_sd_i,
    q025 = pred_q025_i,
    q975 = pred_q975_i
  ) |>
    left_join(
      year_lookup,
      by = "year_idx"
    )
  
  
  ## ---------------------------------------------------------------------------
  ## Diagnostic: confirm the final coordinate grid includes Shannon
  ## ---------------------------------------------------------------------------
  
  pred_const_sf <- st_as_sf(
    df_pred[[const]] |> distinct(X, Y),
    coords = c("X", "Y"),
    crs = kmproj,
    remove = FALSE
  )
  
  shannon_km <- st_transform(shannon, kmproj)
  n_const_shannon <- sum(lengths(st_intersects(pred_const_sf, shannon_km)) > 0)
  
  cat(const, "- prediction locations intersecting Shannon:", n_const_shannon, "\n")
  
  rm(
    pred_model_i,
    fitted_summary_i,
    pred_mean_i,
    pred_sd_i,
    pred_q025_i,
    pred_q975_i,
    pred_const_sf,
    stk.pred.i,
    stack.full.i,
    idx.pred.i
  )
  
  gc(
    verbose = FALSE
  )
}


## =============================================================================
## 26. COMBINE PREDICTIONS
## =============================================================================

prediction_results <- imap_dfr(
  df_pred,
  \(x, constituent) {
    
    if (is.null(x)) {
      return(NULL)
    }
    
    x |>
      mutate(
        constituent = constituent
      )
  }
)

if (nrow(prediction_results) == 0) {
  stop(
    "No spatio-temporal predictions were successfully extracted."
  )
}


## =============================================================================
## 27. CREATE GIS GEOMETRY IN ITM
## =============================================================================

## Convert model coordinates from kilometres to metres
prediction_results_m <- prediction_results |>
  mutate(X_m = X * 1000, Y_m = Y * 1000)

## Create geometry using standard UTM zone 30 (metres)
prediction_gis <- prediction_results_m |>
  st_as_sf(coords = c("X_m", "Y_m"), crs = 32630) |>
  st_transform(2157)

## Extract ITM coordinates
prediction_itm <- st_coordinates(prediction_gis)

prediction_gis <- prediction_gis |>
  mutate(Easting = prediction_itm[, "X"], Northing = prediction_itm[, "Y"])

stopifnot(st_crs(prediction_gis)$epsg == 2157)

## =============================================================================
## 28. CREATE ONE GIS PREDICTION LAYER PER TIDAL FEATURE
## =============================================================================

prediction_layers <- split(
  prediction_gis,
  prediction_gis$constituent
)

prediction_layers <- map(
  prediction_layers,
  \(x) {
    x |>
      select(year_idx, year_mid, Easting, Northing, mean, sd, q025, q975, geometry) |>
      arrange(year_mid, Northing, Easting)
  }
)


## =============================================================================
## 29. FINAL SHANNON COVERAGE CHECK IN ITM GIS OUTPUT
## =============================================================================

shannon_itm <- st_transform(shannon, itmproj)

for (const in names(prediction_layers)) {
  x <- prediction_layers[[const]]
  latest_year <- max(x$year_mid, na.rm = TRUE)
  x_latest <- x |> filter(year_mid == latest_year)
  
  n_final_shannon <- sum(lengths(st_intersects(x_latest, shannon_itm)) > 0)
  cat(const, "- Shannon locations in latest-year GIS layer:", n_final_shannon, "\n")
}


## =============================================================================
## 30. WRITE PREDICTION LAYERS TO GEOPACKAGE
## =============================================================================

for (const in names(prediction_layers)) {
  
  layer_name <- paste0(
    "pred_",
    const
  )
  
  cat(
    "Writing layer:",
    layer_name,
    "\n"
  )
  
  st_write(
    prediction_layers[[const]],
    gpkg_file,
    layer = layer_name,
    append = FALSE,
    quiet = TRUE
  )
}


## =============================================================================
## 31. FINAL GEOPACKAGE CHECKS
## =============================================================================

cat(
  "\n",
  "============================================================\n",
  "GEOPACKAGE EXPORT COMPLETE\n",
  "============================================================\n",
  "\n",
  "File:\n",
  gpkg_file,
  "\n\n",
  sep = ""
)

print(
  st_layers(
    gpkg_file
  )
)

cat(
  "\nStation locations:",
  nrow(stations_gis),
  "\n"
)

cat(
  "Station/year tidal-feature rows:",
  nrow(tidal_features_gis),
  "\n"
)

cat("Prediction locations per year:", n.pred, "
")

cat("Prediction years:", n_years, "\n\n")