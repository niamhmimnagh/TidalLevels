
## =============================================================================
## RASTERISE MI/OPW SPATIO-TEMPORAL PREDICTIONS
## =============================================================================

library(sf)
library(terra)
library(dplyr)
library(here)

## Settings
gpkg_file <- here::here("outputs", "gis", "MI_OPW_tidal_analysis_v2.gpkg")
out_dir <- here::here("outputs", "gis", "rasters")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

const_levels <- c("HAT", "MHWS", "MHWN", "MSL", "MLWN", "MLWS", "LAT")
statistics <- c("mean", "sd", "q025", "q975")

## Original prediction grid resolution (metres)
grid_resolution <- 2000

## Check available layers
available_layers <- sf::st_layers(gpkg_file)$name

## =============================================================================
## RASTERISE EACH CONSTITUENT
## =============================================================================

for (const in const_levels) {
  
  layer_name <- paste0("pred_", const)
  
  if (!layer_name %in% available_layers) {
    warning("Missing layer: ", layer_name)
    next
  }
  
  cat("\nProcessing:", const, "\n")
  
  ## Read prediction layer
  pred <- st_read(gpkg_file, layer = layer_name, quiet = TRUE)
  
  ## Convert ITM coordinates back to original model CRS (UTM 30N)
  pred <- st_transform(pred, 32630)
  
  ## Extract coordinates
  xy <- st_coordinates(pred)
  
  pred <- pred |>
    st_drop_geometry() |>
    mutate(X = xy[, 1], Y = xy[, 2])
  
  years <- sort(unique(pred$year_mid))
  
  ## Create raster template aligned to original 2 km grid
  ## Cell centres are derived from the prediction coordinates.
  x_centres <- sort(unique(round(pred$X, 4)))
  y_centres <- sort(unique(round(pred$Y, 4)))
  
  r_template <- rast(
    xmin = min(x_centres) - grid_resolution / 2,
    xmax = max(x_centres) + grid_resolution / 2,
    ymin = min(y_centres) - grid_resolution / 2,
    ymax = max(y_centres) + grid_resolution / 2,
    resolution = grid_resolution,
    crs = "EPSG:32630"
  )
  
  ## Check original grid alignment
  cells <- cellFromXY(r_template, as.matrix(pred[, c("X", "Y")]))
  xy_check <- xyFromCell(r_template, cells)
  
  if (any(abs(xy_check[, 1] - pred$X) > 1e-3 |
          abs(xy_check[, 2] - pred$Y) > 1e-3)) {
    stop("Prediction points do not align with the 2 km grid for ", const)
  }
  
  ## Rasterise each prediction statistic
  for (stat in statistics) {
    
    cat("  Rasterising:", stat, "\n")
    
    raster_list <- vector("list", length(years))
    
    for (i in seq_along(years)) {
      
      yr <- years[i]
      
      pred_year <- pred |>
        filter(year_mid == yr) |>
        select(X, Y, all_of(stat)) |>
        filter(is.finite(.data[[stat]]))
      
      ## Create raster from the original grid cell centres
      r <- r_template
      
      if (nrow(pred_year) > 0) {
        cell_ids <- cellFromXY(r, as.matrix(pred_year[, c("X", "Y")]))
        
        if (anyDuplicated(cell_ids)) {
          stop("Multiple predictions assigned to one cell: ", const, " ", yr)
        }
        
        r[cell_ids] <- pred_year[[stat]]
      }
      
      names(r) <- paste0("Y", yr)
      raster_list[[i]] <- r
    }
    
    ## Combine all years into a multiband raster
    r_multi <- rast(raster_list)
    names(r_multi) <- paste0("Y", years)
    
    ## Reproject complete raster to Irish Transverse Mercator
    r_itm <- project(r_multi, "EPSG:2157", method = "bilinear")
    names(r_itm) <- names(r_multi)
    
    ## Export GeoTIFF
    out_file <- file.path(out_dir, paste0(const, "_", stat, "_national_2km.tif"))
    
    writeRaster(
      r_itm,
      out_file,
      overwrite = TRUE,
      gdal = c("COMPRESS=DEFLATE", "TILED=YES"),
      datatype = "FLT4S"
    )
    
    cat("  Saved:", basename(out_file), "\n")
    
    rm(r_multi, r_itm, raster_list)
    gc(verbose = FALSE)
  }
}

cat("\nAll tidal prediction rasters exported successfully.\n")
