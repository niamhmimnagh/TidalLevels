library(arrow)
library(tidyverse)
library(patchwork)
theme_set(theme_bw())

fl <- "IrishNationalTideGaugeNetwork_712b_5843_9069"
## As the file has a second row with units, I'll have to skip it, but need the
## col names only first
## col_names <- names(readr::read_csv(paste0(fl, ".csv"), n_max = 0))
## da <- read_csv_arrow(paste0(fl, ".csv"), as_data_frame = FALSE,
##     skip = 2, col_names = col_names)
## write_parquet(da, paste0(fl, ".parquet"))
## rm(da)

da <- read_parquet(paste0(fl, ".parquet"), as_data_frame = FALSE)
class(da)
names(da)

## Check QC flags
da |>
    count(QC_Flag) |>
    collect()

da |>
    count(station_id) |>
    arrange(desc(n)) |>
    collect() |>
    print(n = Inf)

##------------------------------------------------------------------------------
## Select Dublin Port only

da |>
    filter(station_id == "Dublin Port") |>
    distinct(latitude, longitude, station_id, datasourceid) |>
    collect()

db <- da |>
    select(-c(altitude, latitude, longitude, datasourceid)) |>
    filter(station_id == "Dublin Port") |>
    collect()
str(db)
summary(db)

db |>
    count(QC_Flag) |>
    mutate(prop = n/sum(n))

## Make sure it's arranged by time
db <- db |>
    arrange(time)

db <- db |>
    mutate(year = year(time), .before = "station_id")

db |>
    count(year) |>
    arrange(year)

ggplot(db, aes(x = Water_Level_LAT, y = Water_Level_OD_Malin)) +
    geom_point() +
    facet_wrap(~ QC_Flag, scales = "free")
## Clearly QC has to be considered

p1 <- db |>
    filter(QC_Flag == 1) |>
    ggplot() +
    geom_histogram(aes(x = Water_Level_LAT))
p2 <- db |>
    filter(QC_Flag == 1) |>
    ggplot() +
    geom_histogram(aes(x = Water_Level_OD_Malin))
p1 + p2
## Only the scale changes

## TIME
db <- db |>
    mutate(dtime = as.numeric(difftime(time, lag(time), units = "mins")))
sort(unique(db$dtime))

db |>
    filter(time >= as.POSIXct("2007-03-20", tz = "UTC"),
        time < as.POSIXct("2007-03-25", tz = "UTC")) |>
    ggplot(aes(x = time, y = Water_Level_OD_Malin)) +
    geom_line()
## twice-daily (semi-diurnal)

db <- db |>
    mutate(ddays = as.numeric(difftime(time, lag(time), units = "days")))
sort(sort(unique(db$ddays)))

## Try to find the gaps
gaps <- db |>
    mutate(prev_time = lag(time),
        gap_min = as.numeric(difftime(time, prev_time, units = "mins"))) |>
    filter(!is.na(gap_min), gap_min != 5) |>
    select(time, prev_time, gap_min, QC_Flag)

gaps |>
    filter(QC_Flag == 1) |>
    arrange(desc(gap_min)) |>
    head(20)

gaps_ctx <- db |>
    mutate(prev_time = lag(time),
        next_time = lead(time),
        gap_min   = as.numeric(difftime(time, prev_time, units = "mins")),
        next_gap  = as.numeric(difftime(next_time, time, units = "mins"))) |>
    filter(!is.na(gap_min), gap_min != 5) |>
    select(prev_time, time, next_time, gap_min, next_gap, QC_Flag)

gaps_ctx |> arrange(desc(gap_min)) |> head(20)

gap_summary <- gaps |>
    summarise(n_gaps = n(),
        max_gap_min = max(gap_min)) |>
    arrange(desc(n_gaps), desc(max_gap_min))
gap_summary

gaps |>
    mutate(bin = cut(gap_min, breaks = c(0,10,60,180,1440, Inf),
        labels = c("6–10","11–60","61–180","181–1440",">1440"))) |>
    count(bin) |>
    arrange(desc(n))
.Last.value$n |> sum()

db |>
    filter(!is.na(dtime)) |>
    ggplot(aes(x = dtime)) +
    geom_histogram(binwidth = 1000) +
    labs(x = "Interval (minutes)", y = "Count")

## st <- "Dublin Port"
p1 <- db |>
    ## filter(station_id == st) |>
    mutate(is_gap = !is.na(dtime) & dtime != 5) |>
    ggplot(aes(time, Water_Level_OD_Malin)) +
    geom_line() +
    geom_point(data = ~ subset(.x, is_gap), aes(time, Water_Level_OD_Malin),
        color = 2)
p2 <- db |>
    filter(QC_Flag == 1) |>
    mutate(is_gap = !is.na(dtime) & dtime != 5) |>
    ggplot(aes(time, Water_Level_OD_Malin)) +
    geom_line() +
    geom_point(data = ~ subset(.x, is_gap), aes(time, Water_Level_OD_Malin),
        color = 2)
p1/p2

db |>
    mutate(is_gap = !is.na(dtime) & dtime != 5) |>
    ggplot(aes(time, Water_Level_OD_Malin)) +
    geom_line() +
    geom_point(data = ~ subset(.x, is_gap), aes(time, Water_Level_OD_Malin),
        color = 2) +
    facet_wrap(~ QC_Flag, nrow = 2) +
    scale_x_datetime(
        date_breaks = "2 years",
        date_minor_breaks = "1 year",
        date_labels = "%Y")

# compute rounded limits
(lims <- range(db$time, na.rm = TRUE))
(lims <- c(floor_date(lims[1], "2 years"), ceiling_date(lims[2], "2 years")))

db |>
    mutate(is_gap = !is.na(dtime) & dtime != 5) |>
    ggplot(aes(time, Water_Level_LAT)) +
    geom_line() +
    geom_point(
        data = ~ subset(.x, is_gap),
        aes(time, Water_Level_LAT),
        color = 2
    ) +
    facet_wrap(~ QC_Flag, nrow = 2) +
    scale_x_datetime(
        limits = lims,
        date_breaks = "2 years",
        date_labels = "%Y",
        expand = c(0, 0)
    )

# compute rounded limits
(lims <- range(db$time, na.rm = TRUE))
(lims <- c(floor_date(lims[1], "1 year"), ceiling_date(lims[2], "1 year")))

db |>
    mutate(is_gap = !is.na(dtime) & dtime != 5) |>
    ggplot(aes(time, Water_Level_LAT)) +
    geom_line() +
    geom_point(
        data = ~ subset(.x, is_gap),
        aes(time, Water_Level_LAT),
        color = 2
    ) +
    facet_wrap(~ QC_Flag, nrow = 2) +
    scale_x_datetime(
        limits = lims,
        date_breaks = "1 year",
        date_labels = "%Y"
    )

db |>
    count(year, QC_Flag)

db |>
    filter(QC_Flag == 1) |>
    count(year)
## QC_Flag == 1 ends in 2019!!!

##==============================================================================
## Compare with GESLA!
library(geslaR)

## ~15 min
fl <- "gesla-dublin-all"
if(!file.exists(paste0(fl, ".parquet"))) {
    dg <- query_gesla(country = "IRL", site_name = c("DublinPort", "Dublin_Port"),
        use_flag = c(0, 1))
    ## Save on disk
    write_gesla(dg, file_name = fl)
} else {
    dg <- read_gesla(paste0(fl, ".parquet"))
}

names(dg)

dg |>
    group_by(site_name) |>
    count() |>
    collect()

dg |>
    group_by(site_name, year) |>
    count() |>
    arrange(site_name, year) |>
    collect() |>
    print(n = Inf)

dg |>
    group_by(site_name) |>
    count(qc_flag) |>
    collect()

## Remove qc_flag = 3
dg <- dg |>
    filter(qc_flag != 3) |>
    compute()

dg |>
    group_by(site_name) |>
    count(qc_flag) |>
    collect()

dg |>
    group_by(site_name, year) |>
    count() |>
    arrange(site_name, year) |>
    collect() |>
    print(n = Inf)

dg |>
    group_by(site_name) |>
    summarise(
        dmin = min(date_time),
        dmax = max(date_time)
    ) |>
    collect()

dg |>
    group_by(site_name, year, qc_flag) |>
    summarise(m = mean(sea_level, na.rm = TRUE)) |>
    collect() |>
    ggplot(aes(x = year, y = m, color = site_name)) +
    geom_line() +
    facet_wrap(~ qc_flag)

db |>
    group_by(station_id) |>
    summarise(
        dmin = min(time),
        dmax = max(time)
    ) |>
    collect()

db |>
    group_by(year, QC_Flag) |>
    summarise(m = mean(Water_Level_OD_Malin)) |>
    collect() |>
    ggplot(aes(x = year, y = m)) +
    geom_line() +
    facet_wrap(~ QC_Flag)


db |>
    group_by(year, QC_Flag) |>
    count() |>
    arrange(QC_Flag, year)

dg |>
    group_by(year, qc_flag, site_name) |>
    count() |>
    arrange(site_name, year, qc_flag) |>
    collect() |>
    print(n = Inf)
