library(arrow)
library(dplyr)
library(here)

da <- open_dataset(
  here::here(
    "data",
    "merged_tide_gauge_data",
    "iteration_a_qc1"
  ),
  format = "parquet"
)

station_coverage <- da |>
  group_by(st_id) |>
  summarise(
    first_obs = min(time, na.rm = TRUE),
    last_obs = max(time, na.rm = TRUE),
    .groups = "drop"
  ) |>
  collect() |>
  mutate(
    record_years = as.numeric(
      difftime(
        last_obs,
        first_obs,
        units = "days"
      )
    ) / 365.25
  ) |>
  arrange(record_years)

## Stations to exclude
station_coverage |>
  filter(record_years < 2.5)

## Stations slightly short of 3 years but still retained
station_coverage |>
  filter(
    record_years >= 2.5,
    record_years < 3
  )

## All stations
station_coverage