library(tidyverse)

fl_templ <- "tide-stats-template"

run <- tibble(
    station = c("Galway Port"),
    ## station = c("Galway Port", "Dublin Port", "Wexford Harbour",
    ##     "Sligo", "Dingle Harbour"),
    stn = gsub(" ", "", station),
    output_file = sprintf("%s_%s.html", fl_templ, stn),
    execute_params = map(station, \(x) list(station = x))
) |>
    select(-c(station, stn))
run

pwalk(
    .l = run,                      # Dataframe to map over
    .f = quarto::quarto_render,   # Quarto render function
    input = sprintf("%s.qmd", fl_templ),       # Named arguments of .f
    .progress = TRUE              # Optionally, show a progress bar
)


##==============================================================================

library(tidyverse)

fl_templ <- "tide-stats-template-Galway"

run <- tibble(
    station = c("Galway Port"),
    dtime = c(5, 6),
    ## station = c("Galway Port", "Dublin Port", "Wexford Harbour",
    ##     "Sligo", "Dingle Harbour"),
    stn = gsub(" ", "", station),
    ## output_file = sprintf("%s_%s.html", fl_templ, stn),
    output_file = paste0(fl_templ, "_", stn, "_", dtime, "_min.html"),
    execute_params = map2( # Named list of parameters
        station, dtime,
        \(station, dtime) list(station = station, dtime = dtime)
    )
) |>
    select(-c(station, stn))
str(run)

pwalk(
    .l = run |> select(output_file, execute_params),
    .f = quarto::quarto_render,
    input = sprintf("%s.qmd", fl_templ),
    .progress = TRUE
)


##==============================================================================
library(tidyverse)

fl_templ <- "tide-features"

run <- tibble(
    ## station = c("Galway Port"),
    dtime = 10,
    station = c("Galway Port", "Dublin Port", "Wexford Harbour",
        "Sligo", "Dingle Harbour"),
    stn = gsub(" ", "", station),
    ## output_file = sprintf("%s_%s.html", fl_templ, stn),
    output_file = paste0(fl_templ, "_", stn, ".html"),
    execute_params = map2( # Named list of parameters
        station, dtime,
        \(station, dtime) list(station = station, dtime = dtime)
    )
) |>
    select(-c(station, stn))
str(run)

pwalk(
    .l = run |> select(output_file, execute_params),
    .f = quarto::quarto_render,
    input = sprintf("%s.qmd", fl_templ),
    .progress = TRUE
)
