##------------------------------------------------------------------------------
fit_th <- function(x, dto, pred = TRUE, dtime_pred = NULL,
                   t1 = NULL, t2 = NULL, ...) {
    m <- ftide(x = x, dto = dto, ...)
    ## Keep only strictly needed components for prediction (otherwise the object
    ## size is too big)
    m <- list(
        coefficients = m$coefficients,
        call = m$call,
        features1 = m$features1,
        features2 = m$features2,
        msl = m$msl,
        cfmat = m$cfmat,
        apmat = m$apmat,
        fval = m$fval,
        nodal = m$nodal,
        astlon = m$astlon,
        orig = m$orig,
        vn0 = m$vn0
    )
    class(m) <- c("tide", "lm")
    ## Predict for the whole dataset (completed with NAs)
    # The by argument requires time in hours, so 1/60*x should give x minutes
    if(pred) {
        m$pred <- predict(m, from = t1, to = t2, by = 1/60 * dtime_pred)
    }
    return(m)
}

##------------------------------------------------------------------------------
seq_win <- function(x, k, step = k, complete_only = TRUE) {
    stopifnot(is.numeric(x), length(k) == 1L, k >= 1L)
    x <- sort(unique(as.integer(x)))
    if (k > length(x)) return(list())
    year_first <- min(x)
    year_last  <- max(x)
    starts <- seq(from = year_last, to = year_first, by = -step)
    res <- lapply(starts, function(y) seq(y, length.out = k, by = -1))
    if (complete_only) {
        res <- Filter(function(s) all(s %in% x), res)
    }
    res
}

##------------------------------------------------------------------------------
##' @title General summary of a model diagnostics
##'
##' @description This function returns a data frame with the main diagnostics of
##'     an INLA model (DIC, WAIC, lCPO, etc).
##'
##' @details NEEDS UPDATE. The function returns a data frame with the following
##'     columns: DIC, NP.DIC, WAIC, NP.WAIC, CPO, FAIL, NFAIL, CPU. The first
##'     four columns are the main diagnostics of the model, the next two are the
##'     cross-validation diagnostics, and the last one is the CPU time used by
##'     the model.
##'
##' @param mod An object of class \code{inla}
##' @param y The (observed) response variable
##' @param idx The index from INLA stack
##' @param fun Function used to summarise the values
##' @param na.rm Remove NAs (defaults to \code{TRUE})
##' @param KIS Keep It Simple: keep extrictly the main measures only (don't show
##'     number of parameters, fails and CPU time)
##'
##' @return A named vector with the main diagnostics of the model.
##' @author Fernando Mayer
##'
##' @export
mod_summary <- function(mod, y, idx, fun = sum, na.rm = TRUE, KIS = TRUE) {
  res <- c(
    DIC = mod$dic$dic,
    NP.DIC = mod$dic$p.eff,
    WAIC = mod$waic$waic,
    NP.WAIC = mod$waic$p.eff,
    lCPO = -fun(log(mod$cpo$cpo), na.rm = na.rm),
    FAIL = fun(mod$cpo$failure, na.rm = na.rm),
    NFAIL = fun(mod$cpo$failure > 0, na.rm = na.rm),
    lPO = -fun(log(mod$po$po), na.rm = na.rm),
    MSE = mean((y[idx] - mod$summary.fitted.value$mean[idx])^2, na.rm = na.rm),
    MAE = mean(abs(y[idx] - mod$summary.fitted.value$mean[idx]), na.rm = na.rm),
    CPU = unname(mod$cpu.used["Total"])
  )
  
  if (KIS) {
    res <- res[names(res) %in% c("DIC", "WAIC", "lCPO", "lPO", "MSE", "MAE")]
  }
  
  return(res)
}

##------------------------------------------------------------------------------
##' @title Marginal posterior distributions of the parameters
##'
##' @description This function extracts the marginal posterior
##' distributions of the parameters from an INLA model object.
##'
##' @details The function extracts the marginal posterior distributions
##' of the parameters from an INLA model object. The function returns a
##' list of tibbles, one for each parameter. Each tibble contains the
##' marginal posterior distribution of the parameter.
##'
##' @param mod INLA model object
##'
##' @return A list of tibbles, one for each parameter.
##' @author Fernando Mayer
##'
##' @import INLA tidyverse
##' @importFrom purrr map
##' @export
get_marginals <- function(mod) {
    n_fixed <- length(mod$marginals.fixed)
    names_fixed <- names(mod$marginals.fixed)
    n_hyperpar <- length(mod$marginals.hyperpar)
    names_hyperpar <- names(mod$marginals.hyperpar)
    # Create an empty list for results
    res_marg <- vector("list", length = n_fixed + n_hyperpar)
    names(res_marg) <- c(names_fixed, names_hyperpar)
    # Fixed effects (Identity transformation)
    res_marg[names_fixed] <- lapply(mod$marginals.fixed,
        function(m) inla.tmarginal(identity, m))
    # Define transformation functions for hyperparameters
    hyperpar_transforms <- list(
        "Precision for the Gaussian observations" = \(x) 1/sqrt(x),
        "Precision for id2" = \(x) 1/sqrt(x),
        "Precision for id3" = \(x) 1/sqrt(x),
        "Rho for id3" = \(x) x
    )
    # Apply transformations dynamically
    for (nh in names_hyperpar) {
        ## %||% = if(is.null(x)) y else x
        transform_fun <- hyperpar_transforms[[nh]] %||% function(x) x
        res_marg[[nh]] <- inla.tmarginal(transform_fun,
            mod$marginals.hyperpar[[nh]],
            method = "linear")
    }
    # Define renaming mapping
    rename_map <- c(
        "Precision for the Gaussian observations" = "Sigma",
        "Precision for id2" = "Sigma_i",
        "Precision for id3" = "Sigma_t",
        "Rho for id3" = "Rho"
    )
    # Rename list elements based on mapping
    names(res_marg) <- ifelse(
        names(res_marg) %in% names(rename_map),
        rename_map[names(res_marg)],
        names(res_marg)
    )
    # Convert to tibbles
    res_marg <- purrr::map(res_marg, as_tibble)
    # Assign class
    class(res_marg) <- c("res_marg", class(res_marg))
    return(res_marg)
}

##' @title Summary of the marginal posterior distributions of the
##' parameters
##'
##' @description This function summarizes the marginal posterior
##' distributions of the parameters.
##'
##' @details The function summarizes the marginal posterior
##' distributions of the parameters. The function returns a table with
##' the mean, median, mode, standard deviation, 2.5% quantile, and 97.5%
##' quantile of the marginal posterior distributions of the parameters.
##'
##' @param object A list of tibbles, one for each parameter.
##' @param round Number of decimal places to round the results.
##' @param ... Additional arguments to be passed to the summary function.
##'
##' @return A table with the mean, median, mode, standard deviation,
##' 2.5% quantile, and 97.5% quantile of the marginal posterior
##' distributions of the parameters.
##' @author Fernando Mayer
##'
##' @import INLA tidyverse
##' @importFrom purrr map
##' @method summary res_marg
##' @export
summary.res_marg <- function(object, round = NULL, ...) {
    marginals <- object
    ## Summary table
    res_tab <- purrr::map(marginals,
        function(x) inla.zmarginal(x, silent = TRUE))
    res_tab <- purrr::map(res_tab, as_tibble) |>
        bind_rows(.id = "parameter")
    res_tab <- res_tab[, -c(5, 7)]
    if(!is.null(round)) {
        res_tab <- res_tab |>
            mutate(across(where(is.numeric), \(x) round(x, round)))
    }
    ## UseMethod("summary")
    return(res_tab)
}

##' @title Plot marginal posterior distributions of the parameters
##'
##' @description This function plots the marginal posterior
##' distributions of the parameters.
##'
##' @details The function plots the marginal posterior distributions of
##' the parameters. The function returns a plot with the marginal
##' posterior distributions of the parameters.
##'
##' @param x A list of tibbles, one for each parameter.
##' @param ... Additional arguments to be passed to the plot function.
##'
##' @return A plot with the marginal posterior distributions of the
##' parameters.
##'
##' @author Fernando Mayer
##'
##' @import ggplot2 dplyr
##' @importFrom forcats fct_relevel
##' @importFrom purrr map
##' @method plot res_marg
##' @export
plot.res_marg <- function(x, ...) {
    marginals <- x
    res_marg <- purrr::map(marginals, as_tibble) |>
        bind_rows(.id = "Parameter") |>
        mutate(
            Parameter =
                forcats::fct_relevel(Parameter, names(marginals))
        )
    ggplot(res_marg, aes(x = x, y = y)) +
        geom_line() +
        facet_wrap(~ Parameter, scales = "free",
            nrow = nrow(res_marg)) +
        xlab("") + ylab("Density")
    ## UseMethod("plot")
}

##------------------------------------------------------------------------------
##' @title Get fit and residuals from INLA model
##'
##' @description Get fit and residuals from INLA model for a given
##' dataset and model.
##'
##' @details This function extracts the fitted values and residuals from
##' an INLA model for a given dataset. The function can also be used to
##' extract the fitted values and residuals from a given function.
##'
##' @param mod An INLA model object.
##' @param y A numeric vector of observed values.
##' @param idat A numeric vector of indices.
##' @param fun A function to extract the fitted values and residuals.
##' @param std A logical value indicating whether to standardize the
##' residuals. Default is \code{TRUE}.
##'
##' @return A data frame with the fitted values and residuals.
##' @author Fernando Mayer
##'
##' @importFrom INLA inla.tmarginal inla.zmarginal
##' @importFrom parallel mclapply detectCores
##' @export
get_fit_res <- function(mod, y, idat, fun = NULL, std = TRUE){
    ncores <- as.integer(parallel::detectCores() - 2)
    if(is.null(fun)){
        fit.mean <- mod$summary.fitted.values$mean[idat]
        fit.sd <- mod$summary.fitted.values$sd[idat]
        fit.q <- mod$summary.fitted.values[idat, 3:5]
        names(fit.q) <- c("Q1.fit", "Q2.fit", "Q3.fit")
    } else{
        marg <- mod$marginals.fitted.values
        tmarg <- parallel::mclapply(marg, function(x){
            INLA::inla.tmarginal(fun = fun, x)
        }, mc.cores = ncores)
        res.marg <- parallel::mclapply(tmarg, function(x){
            unlist(INLA::inla.zmarginal(x, silent = TRUE)[1:2])
        }, mc.cores = ncores)
        ## res.marg <- t(res.marg)
        res.marg <- do.call(rbind, res.marg)
        res.marg <- as.data.frame(res.marg,
            row.names = 1:nrow(res.marg))
        fit.mean <- res.marg$mean
        fit.sd <- res.marg$sd
    }
    if(std){
        res <- (y - fit.mean)/fit.sd
    } else{
        res <- (y - fit.mean)
    }
    out <- cbind(y, fit.mean, fit.q, res = res)
    rownames(out) <- NULL
    class(out) <- c("fit_res", class(out))
    return(out)
}

##------------------------------------------------------------------------------
##' @title Generate all combinations of INLA model formulas
##'
##' @description This function generates all sensible combinations of INLA model
##'     formulas based on specified components (intercept, linear trend, and
##'     random effects).
##'
##' @param response Character string with the response variable name. Default is
##'     "value".
##' @param linear_trend Character vector with linear trend options. Default is
##'     c("", "id") where "" means no linear trend.
##' @param random_effects Named list of random effect specifications. Each
##'     element should be a character string with the random effect formula or
##'     "" for no random effect.
##' @param intercept Logical. Include intercept? Default is TRUE.
##'
##' @return A named list of formulas.
##' @author Fernando Mayer
##'
##' @export
##'
##' @examples
##' # Default usage
##' forms <- generate_model_formulas()
##'
##' # Custom random effects
##' re <- list(
##'     none = "",
##'     iid = "f(id2, model = 'iid')",
##'     ar1 = "f(id3, model = 'ar1')"
##' )
##' forms <- generate_model_formulas(random_effects = re)
gen_model_formulas <- function(response = "value", linear_trend = c("",
    "id"), random_effects = list(norandom = "", iid = "f(id2, model = 'iid')", ar1
    = "f(id3, model = 'ar1')", rw1 = "f(id3, model = 'rw1')", iid_ar1 = "f(id2,
    model = 'iid') + f(id3, model = 'ar1')", iid_rw1 = "f(id2, model = 'iid') +
    f(id3, model = 'rw1')" ), intercept = TRUE) {
    # Initialize results list
    form_list <- list()
    counter <- 1
    # Intercept component
    int_part <- if (intercept) "1" else "0"
    # Loop through all combinations
    for (trend in linear_trend) {
        for (re_name in names(random_effects)) {
            re <- random_effects[[re_name]]
            # Build formula parts
            parts <- c(int_part, trend, re)
            parts <- parts[parts != ""]  # Remove empty strings
            # Create formula
            formula_str <- paste(response, "~", paste(parts, collapse = " + "))
            form_list[[counter]] <- as.formula(formula_str)
            # Create descriptive name
            trend_label <- ifelse(trend == "", "int", "trend")
            names(form_list)[counter] <- paste0("m", counter, "_",
                trend_label, "_", re_name)
            counter <- counter + 1
        }
    }
    return(form_list)
}


##==============================================================================
##' @title Generate spatial and spatio-temporal INLA model formulas
##'
##' @description This function generates formulas for spatial and spatio-temporal
##'     models using SPDE specifications for barrier models in INLA.
##'
##' @param response Character string with the response variable name. Default is
##'     "y".
##' @param fixed_effects Character vector with fixed effect terms. Default is
##'     c("m", "trendc") for intercept and centered trend.
##' @param spatial_field Character string with the spatial random effect
##'     specification. Default uses the barrier model setup.
##' @param temporal_models Character vector with temporal model types. Options:
##'     "none" (spatial only), "iid", "ar1", "rw1". Default is all.
##' @param w_index Character string with the name of the space-time index.
##'     Default is "w.index".
##' @param intercept Logical. Include intercept (m)? Default is TRUE.
##' @param trend Logical. Include linear trend? Default is TRUE.
##'
##' @return A named list of formulas with spatial models first, then
##'     spatio-temporal models.
##' @author Fernando Mayer
##'
##' @export
##'
##' @examples
##' # Default usage
##' forms <- gen_st_model_formulas()
##'
##' # Only spatial and AR1
##' forms <- gen_st_model_formulas(temporal_models = c("none", "ar1"))
##'
##' # No trend
##' forms <- gen_st_model_formulas(trend = FALSE)
gen_st_model_formulas <- function(
                                  response = "y",
                                  fixed_effects = c("m", "trendc"),
                                  spatial_field = "f(w, model = barrier_model)",
                                  temporal_models = c("none", "iid", "ar1", "rw1"),
                                  w_index = "w.index",
                                  intercept = TRUE,
                                  trend = TRUE
                                  ) {
    # Initialize results list
    form_list <- list()
    counter <- 1
    # Build fixed effects part
    fixed_parts <- c()
    if (intercept && "m" %in% fixed_effects) {
        fixed_parts <- c(fixed_parts, "m")
    }
    if (trend && "trendc" %in% fixed_effects) {
        fixed_parts <- c(fixed_parts, "trendc")
    }
    # If no fixed effects, use 0
    if (length(fixed_parts) == 0) {
        fixed_str <- "0"
    } else {
        fixed_str <- paste(fixed_parts, collapse = " + ")
    }
    # --- SPATIAL MODELS FIRST ---
    # Pure spatial model
    if ("none" %in% temporal_models) {
        formula_str <- paste(response, "~", fixed_str, "+", spatial_field)
        form_list[[counter]] <- as.formula(formula_str)
        names(form_list)[counter] <- paste0("m", counter, "_spatial")
        counter <- counter + 1
    }
    # --- SPATIO-TEMPORAL MODELS ---
    temporal_specs <- list(
        iid = sprintf("f(%s, model = 'iid')", w_index),
        ar1 = sprintf("f(%s, model = 'ar1', constr = FALSE, hyper = prior_rho)", w_index),
        rw1 = sprintf("f(%s, model = 'rw1', constr = FALSE, scale.model = TRUE)", w_index)
    )
    for (temp_name in temporal_models) {
        if (temp_name == "none") next  # Already handled above
        if (!temp_name %in% names(temporal_specs)) {
            warning("Unknown temporal model: ", temp_name, ". Skipping.")
            next
        }
        temp_spec <- temporal_specs[[temp_name]]
        # Build formula: fixed + spatial + temporal
        formula_str <- paste(response, "~", fixed_str, "+",
            spatial_field, "+", temp_spec)
        form_list[[counter]] <- as.formula(formula_str)
        # Create descriptive name
        trend_label <- if (trend) "trend" else "int"
        names(form_list)[counter] <- paste0("m", counter, "_st_",
            trend_label, "_", temp_name)
        counter <- counter + 1
    }
    return(form_list)
}

##==============================================================================
##' @title Generate spatial and spatio-temporal formulas for extreme value models
##'
##' @description This function generates formulas for spatial and spatio-temporal
##'     models with barrier SPDE, including options for intercept-only, trend,
##'     spatial, and spatio-temporal structures. Each model is returned as a list
##'     with the formula and number of hyperparameters.
##'
##' @param response Character string with the response variable name. Default is
##'     "y".
##' @param intercept Character string for intercept term. Default is "m".
##' @param trend Character string for trend term. Default is "trendc".
##' @param spatial_field Character string for spatial field name. Default is "w".
##' @param barrier_model Character string with barrier model object name. Default
##'     is "barrier_model".
##' @param temporal_group Character string for temporal grouping variable. Default
##'     is "w.group".
##' @param temporal_model Character string for temporal model type. Default is
##'     "ar1".
##' @param prior_rho Character string for prior on temporal correlation. Default
##'     is "prior_rho".
##' @param include_baseline Logical. Include baseline models (m0, m1) without
##'     spatial effects? Default is TRUE.
##' @param include_spatial Logical. Include spatial models (m0s, m1s)? Default is
##'     TRUE.
##' @param include_spatiotemporal Logical. Include spatio-temporal models (m0st,
##'     m1st)? Default is TRUE.
##'
##' @return A named list where each element is a list with two components:
##'     \item{form}{The model formula}
##'     \item{nhyper}{Number of hyperparameters (1, 1, 3, 3, 4, 4)}
##' @author Fernando Mayer
##'
##' @export
##'
##' @examples
##' # Default usage - all models
##' forms <- gen_barrier_model_formulas()
##' forms$m0$form  # Access formula
##' forms$m0$nhyper  # Access number of hyperparameters
##'
##' # Only spatial and spatio-temporal
##' forms <- gen_barrier_model_formulas(include_baseline = FALSE)
gen_barrier_model_formulas_DONT <- function(
                                       response = "y",
                                       intercept = "m",
                                       trend = "trendc",
                                       spatial_field = "w",
                                       barrier_model = "barrier_model",
                                       temporal_group = "w.group",
                                       temporal_model = "ar1",
                                       prior_rho = "prior_rho",
                                       include_baseline = TRUE,
                                       include_spatial = TRUE,
                                       include_spatiotemporal = TRUE
                                       ) {
    form_list <- list()
    # --- Baseline models (no spatial effects) ---
    if (include_baseline) {
        # m0: intercept only
        form_list$m0 <- list(
            form = paste(response, "~ -1 +", intercept),
            nhyper = 1
        )
        # m1: intercept + trend
        form_list$m1 <- list(
            form = paste(response, "~ -1 +", intercept, "+", trend),
            nhyper = 1
        )
    }
    # --- Spatial models ---
    if (include_spatial) {
        # m0s: intercept + spatial
        spatial_spec <- sprintf("f(%s, model = %s)", spatial_field, barrier_model)
        form_list$m0s <- list(
            form = paste(response, "~ -1 +", intercept, "+", spatial_spec),
            nhyper = 3
        )
        # m1s: intercept + trend + spatial
        form_list$m1s <- list(
            form = paste(response, "~ -1 +", intercept, "+",
                trend, "+", spatial_spec),
            nhyper = 3
        )
    }
    # --- Spatio-temporal models ---
    if (include_spatiotemporal) {
        # Build spatio-temporal specification
        st_spec <- sprintf(
            "f(%s, model = %s, group = %s, control.group = list(model = '%s', hyper = %s))",
            spatial_field, barrier_model, temporal_group, temporal_model, prior_rho
        )
        # m0st: intercept + spatio-temporal
        form_list$m0st <- list(
            form = paste(response, "~ -1 +", intercept, "+", st_spec),
            nhyper = 4
        )
        # m1st: intercept + trend + spatio-temporal
        form_list$m1st <- list(
            form = paste(response, "~ -1 +", intercept, "+",
                trend, "+", st_spec),
            nhyper = 4
        )
    }
    return(form_list)
}

##------------------------------------------------------------------------------
gen_barrier_model_formulas <- function() {
    fms <- list(
        m0 = list(form = y ~ -1 + m,
            nhyper = 1),
        m1 = list(form = y ~ -1 + m + trendc,
            nhyper = 1),
        m0s = list(form = y ~ -1 + m +
                       f(w, model = barrier_model),
            nhyper = 3),
        m1s = list(form = y ~ -1 + m + trendc +
                       f(w, model = barrier_model),
            nhyper = 3),
        m0st = list(form = y ~ -1 + m +
                        f(w, model = barrier_model,
                            group = w.group,
                            control.group = list(model = 'ar1',
                                hyper = prior_rho)),
            nhyper = 4),
        m1st = list(form = y ~ -1 + m + trendc +
                        f(w, model = barrier_model,
                            group = w.group,
                            control.group = list(model = 'ar1',
                                hyper = prior_rho)),
            nhyper = 4)
    )
    return(fms)
}

##------------------------------------------------------------------------------
## Functions to use on tabsets

mod_summ <- function(mod, data) {
    ## Group CV
    gcv <- map(mod, \(x) inla.group.cv(x, num.level.sets = 2, verbose = TRUE))
    gcv <- purrr::map_dbl(gcv, \(x) - sum(log(x$cv), na.rm = TRUE))
    ## Model summary
    summ <- map(mod, \(x) mod_summary(x, y = data$value, idx = data$id)) |>
        bind_rows(.id = "mod") |>
        mutate(GCV = gcv)
    summ
}

res_marg_bind <- function(mod) {
    ## Get marginals
    res_marg <- map(mod, get_marginals)
    ## Make the right parameter order
    param_levels <- c("(Intercept)", "id", "Sigma", "Sigma_i", "Sigma_t", "Rho")
    ## Arrange for plotting
    res_marg_bind <- bind_rows(res_marg, .id = "model") |>
        pivot_longer(cols = -model,
            names_to = "parameter",
            values_to = "dens") |>
        unnest(dens)
    res_marg_bind <- res_marg_bind |>
        mutate(parameter = fct_relevel(parameter, param_levels),
            model = fct_relevel(model, names(mod)))
}

mod_marg <- function(mod, constituent) {
  
  res_marg_bind <- res_marg_bind(mod)
  
  ## Use ggridges to plot density ridges
  res_marg_bind |>
    mutate(
      model = sub("_.*", "", model),
      model = fct_relevel(
        model,
        paste0("m", seq_along(mod))
      )
    ) |>
    ggplot(
      aes(
        x = x,
        y = model,
        height = y,
        fill = model
      )
    ) +
    ggridges::geom_density_ridges(
      stat = "identity",
      scale = 1,
      alpha = .7
    ) +
    facet_wrap(
      ~ parameter,
      scales = "free"
    ) +
    labs(
      x = "Value",
      y = "Model",
      fill = "Model",
      title = sprintf(
        "Posterior marginals for [%s] %s",
        constituent,
        params$station
      )
    ) +
    scale_fill_viridis_d() +
    theme(
      legend.position = "none"
    )
}

fit_res_df <- function(mod, data) {
    ## Get fitted and residuals
    fit_res <- map(mod, \(x) get_fit_res(x, y = data$value,
        idat = data$id, std = FALSE))
    ## Arrange for plotting
    fit_res_df <- map(fit_res, \(x) {
        as_tibble(x) |>
            mutate(
                id = data$id,
                year_mid = data$year_mid
            )
    }) |>
        bind_rows(.id = "model")
    fit_res_df <- fit_res_df |>
        mutate(model = fct_relevel(model, names(mod)))
    return(fit_res_df)
}

mod_pred <- function(mod, data, constituent) {
  
  fit_res_df <- fit_res_df(
    mod,
    data
  )
  
  ## Plot observed vs fitted over time
  ggplot(
    fit_res_df,
    aes(x = year_mid)
  ) +
    geom_point(
      aes(y = y),
      color = "black",
      alpha = 0.6
    ) +
    geom_line(
      aes(
        y = fit.mean,
        color = model
      )
    ) +
    geom_ribbon(
      aes(
        ymin = Q1.fit,
        ymax = Q3.fit,
        fill = model
      ),
      alpha = 0.2
    ) +
    facet_wrap(
      ~ model,
      ncol = 2,
      scales = "free_y"
    ) +
    labs(
      x = "Year",
      y = "Value",
      title = sprintf(
        "Observed and predicted values [%s] %s",
        constituent,
        params$station
      )
    ) +
    theme(
      legend.position = "none"
    ) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(
        n = 20
      )
    ) +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5
      )
    )
}

mod_obs_fit <- function(mod, data, constituent) {
  
  fit_res_df <- fit_res_df(
    mod,
    data
  )
  
  ## Observed vs fitted
  fit_res_df |>
    ggplot() +
    geom_point(
      aes(
        x = y,
        y = fit.mean
      ),
      color = "black",
      alpha = 0.6
    ) +
    facet_wrap(
      ~ model,
      ncol = 2
    ) +
    labs(
      x = "Observed",
      y = "Fitted",
      title = sprintf(
        "Observed vs Fitted Values [%s] %s",
        constituent,
        params$station
      )
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      alpha = 0.6,
      col = 2
    )
}

mod_obs_res <- function(mod, data, constituent) {
  
  fit_res_df <- fit_res_df(
    mod,
    data
  )
  
  ## Observed vs residuals
  ggplot(
    fit_res_df,
    aes(
      x = y,
      y = res
    )
  ) +
    geom_point(
      alpha = 0.6
    ) +
    facet_wrap(
      ~ model,
      ncol = 2
    ) +
    labs(
      x = "Observed",
      y = "Residuals",
      title = sprintf(
        "Observed vs Residuals [%s] %s",
        constituent,
        params$station
      )
    ) +
    geom_hline(
      yintercept = 0,
      col = 2,
      linetype = "dashed",
      alpha = 0.6
    )
}

out_table <- function(tab) {
  tab |>
    DT::datatable(
      class = "compact",
      rownames = FALSE,
      options = list(
        ordering = TRUE,
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        dom = "t"
      )
    ) |>
    DT::formatRound(
      columns = setdiff(names(tab)[2:ncol(tab)], c("MSE", "MAE")),
      digits = 3
    ) |>
    DT::formatRound(
      columns = intersect(c("MSE", "MAE"), names(tab)),
      digits = 3
    )
}

##==============================================================================
plot_time_range <- function(data,
                            year_start_col = "year_start",
                            year_end_col = "year_end",
                            id_col = "st_id",
                            title = "Time range of observations by station") {
    id_sym <- sym(id_col)
    data |>
        group_by(!!id_sym) |>
        summarise(
            min_year = make_datetime(min(.data[[year_start_col]], na.rm = TRUE)),
            max_year = make_datetime(max(.data[[year_end_col]], na.rm = TRUE)),
            .groups = "drop"
        ) |>
        mutate(!!id_sym := fct_reorder(!!id_sym, min_year)) |>
        ggplot(aes(y = !!id_sym)) +
        geom_segment(aes(x = min_year, xend = max_year, yend = !!id_sym)) +
        geom_point(aes(x = min_year), color = "blue") +
        geom_point(aes(x = max_year), color = "red") +
        scale_x_datetime(
            date_breaks = "2 year",
            date_minor_breaks = "1 year",
            date_labels = "%Y") +
        labs(x = "Year", y = "", title = title)
}

##==============================================================================
plot_feat_time <- function(data, year_col = "year_mid", value_col = "value",
                           color_col = "constituent", facet_col = "st_id",
                           title = "Estimated features over the years",
                           x_label = "Year (middle)", y_label = "Value",
                           color_label = "Constituent", scales = "fixed",
                           ncol = 4) {
  
  data |>
    bind_rows() |>
    filter(!is.na(.data[[value_col]])) |>
    ggplot(aes(x = make_datetime(.data[[year_col]]), y = .data[[value_col]], color = .data[[color_col]])) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    facet_wrap(as.formula(paste("~", facet_col)), scales = scales, ncol = ncol) +
    scale_x_datetime(date_breaks = "4 years", date_minor_breaks = "1 year", date_labels = "%Y") +
    labs(x = x_label, y = y_label, title = title, color = color_label) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 8),
      legend.position = "right"
    )
}

##==============================================================================
##==============================================================================
plot_feat_dist <- function(data,
                           value_col = "value",
                           facet_col = "constituent",
                           bins = 30,
                           title = "Distribution of features",
                           scales = "free") {
  
  data |>
    bind_rows() |>
    filter(!is.na(.data[[value_col]])) |>
    ggplot(aes(x = .data[[value_col]])) +
    geom_histogram(bins = bins) +
    facet_grid(as.formula(paste("~", facet_col)), scales = scales) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    labs(x = "Value", y = "Count", title = title) +
    theme(
      axis.text.x = element_text(size = 8),
      strip.text = element_text(size = 10)
    )
}

##==============================================================================
add_coords_to_sf <- function(sf_data, coord_names = c("lon", "lat")) {
    sf_data |>
        st_coordinates() |>
        as_tibble() |>
        setNames(coord_names) |>
        bind_cols(sf_data) |>
        st_as_sf()
}

##==============================================================================
center_scale_vars <- function(data,
                              trend_col = "year_idx",
                              lat_col = "lat",
                              lon_col = "lon",
                              return_params = FALSE) {
    # Calculate parameters
    params <- list(
        trend_center = mean(data[[trend_col]], na.rm = TRUE),
        trend_sd = sd(data[[trend_col]], na.rm = TRUE),
        lat_center = mean(data[[lat_col]], na.rm = TRUE),
        lon_center = mean(data[[lon_col]], na.rm = TRUE),
        lat_sd = sd(data[[lat_col]], na.rm = TRUE),
        lon_sd = sd(data[[lon_col]], na.rm = TRUE)
    )
    # Apply centering
    result <- tibble(
        trendc = data[[trend_col]] - params$trend_center,
        latc = (data[[lat_col]] - params$lat_center) / params$lat_sd,
        lonc = (data[[lon_col]] - params$lon_center) / params$lon_sd
    )
    if (return_params) {
        return(list(data = result, params = params))
    } else {
        return(result)
    }
}


##==============================================================================
##' @title Create INLA stacks for spatial and spatio-temporal models
##'
##' @description Creates spatial and spatio-temporal INLA stacks for a single
##'     response variable.
##'
##' @param data Data frame with the response variable.
##' @param spatial_data Spatial data frame with coordinates and covariates.
##' @param response Character string with response variable name. Default is
##'     "value".
##' @param trend Character string with trend variable name. Default is "trendc".
##' @param barrier_model Barrier model object.
##' @param A_spatial Projector matrix for spatial model.
##' @param A_spatiotemporal Projector matrix for spatio-temporal model.
##' @param w_index Space-time index object.
##'
##' @return A list with two elements: spatial and spatiotemporal stacks.
##' @author Fernando Mayer
##'
##' @export
create_inla_stacks <- function(
                               data,
                               spatial_data,
                               response = "value",
                               trend = "trendc",
                               barrier_model,
                               A_spatial,
                               A_spatiotemporal,
                               w_index
                               ) {
    # Extract response and trend
    y_data <- data[[response]]
    trend_data <- spatial_data[[trend]]
    # Spatial stack
    stk_s <- inla.stack(
        data = list(y = y_data),
        effects = list(
            w = 1:barrier_model$f$n,
            data.frame(
                m = 1,
                trendc = trend_data
            )
        ),
        A = list(A_spatial, 1),
        tag = "est.s"
    )
    # Spatio-temporal stack
    stk_st <- inla.stack(
        data = list(y = y_data),
        effects = list(
            w_index,
            data.frame(
                m = 1,
                trendc = trend_data
            )
        ),
        A = list(A_spatiotemporal, 1),
        tag = "est.st"
    )
    return(list(
        spatial = stk_s,
        spatiotemporal = stk_st
    ))
}

##==============================================================================
inla_stack_s <- function(y, barrier, trend, A, tag) {
    stk <- inla.stack(
        data = list(y = y),
        effects = list(
            w = 1:barrier$f$n,
            data.frame(
                m = 1,
                trendc = trend
            )
        ),
        A = list(A, 1),
        tag = tag
    )
}

inla_stack_st <- function(y, w, trend, A, tag) {
    stk <- inla.stack(
        data = list(y = y),
        effects = list(
            w,
            data.frame(
                m = 1,
                trendc = trend
            )
        ),
        A = list(A, 1),
        tag = tag
    )
}

##==============================================================================
inla_st <- function(formula, data, stk, st = TRUE, ...) {
    inla_family <- "gaussian"
    inla_control_compute <- list(cpo = TRUE, waic = TRUE, dic = TRUE,
        config = TRUE)
    inla_verbose <- TRUE
    nt <- ifelse(st, formula$nhyper + 1, 2 * formula$nhyper)
    inla_control_pred <- list(A = inla.stack.A(stk), compute = TRUE, link = 1)
    m <- inla(
        formula = formula$form,
        data = data,
        family = inla_family,
        control.compute = inla_control_compute,
        control.predictor = inla_control_pred,
        verbose = inla_verbose,
        num.threads = nt,
        ...
    )
    return(m)
}

##==============================================================================
mod_summ_st <- function(mod, data, idat) {
    ## Group CV
    ## gcv <- map(mod, \(x) inla.group.cv(x, num.level.sets = sets, verbose = FALSE))
    ## gcv <- purrr::map_dbl(gcv, \(x) - sum(log(x$cv), na.rm = TRUE))
    ## Model summary
    summ <- map(mod, \(x) mod_summary(x, y = data$value, idx = idat)) |>
        bind_rows(.id = "mod")
    summ
}

##------------------------------------------------------------------------------
get_marginals_st <- function(mod) {
    n_fixed <- length(mod$marginals.fixed)
    names_fixed <- names(mod$marginals.fixed)
    n_hyperpar <- length(mod$marginals.hyperpar)
    names_hyperpar <- names(mod$marginals.hyperpar)
    # Create an empty list for results
    res_marg <- vector("list", length = n_fixed + n_hyperpar)
    names(res_marg) <- c(names_fixed, names_hyperpar)
    # Fixed effects (Identity transformation)
    res_marg[names_fixed] <- lapply(mod$marginals.fixed,
        function(m) inla.tmarginal(identity, m))
    # Define transformation functions for hyperparameters
    hyperpar_transforms <- list(
        "Precision for the Gaussian observations" = \(x) 1/sqrt(x),
        "Theta1 for w" = \(x) exp(x),
        "Theta2 for w" = \(x) exp(x),
        "GroupRho for w" = \(x) x
    )
    # Apply transformations dynamically
    for (nh in names_hyperpar) {
        ## %||% = if(is.null(x)) y else x
        transform_fun <- hyperpar_transforms[[nh]] %||% function(x) x
        res_marg[[nh]] <- inla.tmarginal(transform_fun,
            mod$marginals.hyperpar[[nh]],
            method = "linear")
    }
    # Define renaming mapping
    rename_map <- c(
        "Precision for the Gaussian observations" = "Sigma",
        "Theta1 for w" = "rangeM",
        "Theta2 for w" = "sigmaM",
        "GroupRho for w" = "rho"
    )
    # Rename list elements based on mapping
    names(res_marg) <- ifelse(
        names(res_marg) %in% names(rename_map),
        rename_map[names(res_marg)],
        names(res_marg)
    )
    # Convert to tibbles
    res_marg <- purrr::map(res_marg, as_tibble)
    # Assign class
    class(res_marg) <- c("res_marg", class(res_marg))
    return(res_marg)
}

res_marg_bind_st <- function(mod) {
    ## Get marginals
    res_marg <- map(mod, get_marginals_st)
    ## Make the right parameter order
    param_levels <- c("m", "trendc", "Sigma", "rangeM", "sigmaM", "rho")
    ## Arrange for plotting
    res_marg_bind <- bind_rows(res_marg, .id = "model") |>
        pivot_longer(cols = -model,
            names_to = "parameter",
            values_to = "dens") |>
        unnest(dens)
    res_marg_bind <- res_marg_bind |>
        mutate(parameter = fct_relevel(parameter, param_levels),
            model = fct_relevel(model, names(mod)))
    return(res_marg_bind)
}

mod_marg_st <- function(mod) {
    res_marg_bind <- res_marg_bind_st(mod)
    mlevels <- names(mod)
    ## Use ggridges to plot density ridges
    res_marg_bind |>
        mutate(model = fct_relevel(model, mlevels)) |>
        ggplot(aes(x = x, y = model, height = y, fill = model)) +
        ggridges::geom_density_ridges(stat = "identity", scale = 1, alpha = .7) +
        facet_wrap(~ parameter, scales = "free") +
        labs(x = "Value", y = "Model", fill = "Model",
            title = sprintf("Posterior marginals for [%s]", const)) +
        scale_fill_viridis_d() +
        theme(legend.position = "none")
}

get_fit_res_st <- function(mod, y, idat, fun = NULL, std = TRUE){
    ncores <- as.integer(parallel::detectCores() - 2)
    if(is.null(fun)){
        fit.mean <- mod$summary.fitted.values$mean[idat]
        fit.sd <- mod$summary.fitted.values$sd[idat]
        fit.q <- mod$summary.fitted.values[idat, 3:5]
        names(fit.q) <- c("Q1.fit", "Q2.fit", "Q3.fit")
    } else{
        marg <- mod$marginals.fitted.values
        tmarg <- parallel::mclapply(marg, function(x){
            INLA::inla.tmarginal(fun = fun, x)
        }, mc.cores = ncores)
        res.marg <- parallel::mclapply(tmarg, function(x){
            unlist(INLA::inla.zmarginal(x, silent = TRUE)[1:2])
        }, mc.cores = ncores)
        ## res.marg <- t(res.marg)
        res.marg <- do.call(rbind, res.marg)
        res.marg <- as.data.frame(res.marg,
            row.names = 1:nrow(res.marg))
        fit.mean <- res.marg$mean
        fit.sd <- res.marg$sd
    }
    if(std){
        res <- (y - fit.mean)/fit.sd
    } else{
        res <- (y - fit.mean)
    }
    out <- cbind(y, fit.mean, fit.q, res = res)
    rownames(out) <- NULL
    class(out) <- c("fit_res", class(out))
    return(out)
}


fit_res_df_st <- function(mod, data, idat, ...) {
  
  ## Ensure every model has a name
  if (is.null(names(mod)) || any(names(mod) == "")) {
    names(mod) <- paste0("model_",seq_along(mod) )
  }
  
  ## Extract fitted values separately for each model
  fit_res <- purrr::imap(mod,function(model_object, model_name) {
    fitted_summary <- model_object$summary.fitted.values[ idat,,drop = FALSE]
    
    tibble::tibble(
      model = model_name,
      fit.mean = fitted_summary$mean,
      Q1.fit = fitted_summary$`0.025quant`,
      Q3.fit = fitted_summary$`0.975quant`
    )
  }
  )
  
  ## Combine all models and retain the model identifier
  fit_res_df <- dplyr::bind_rows(fit_res)
  
  ## Repeat the observation data once for each fitted model
  data_for_plot <- data |>dplyr::select(dplyr::any_of(c("year_mid","st_id","value")))
  
  fit_res_df <- fit_res_df |>
    dplyr::group_split(model,.keep = TRUE) |>
    purrr::map_dfr(function(model_fit) {
      if (nrow(model_fit) != nrow(data_for_plot)) {
        stop(
          "The number of fitted values does not match the ",
          "number of rows in the constituent data."
        )
      }
      
      dplyr::bind_cols(data_for_plot, model_fit)
    }
    ) |>
    dplyr::rename(y = value) |>
    dplyr::mutate(model = factor(model,levels = names(mod)))
  
  fit_res_df
}


  
  mod_pred_st <- function(mod,data,idat,constituent = NULL,...) {
    
    if (is.null(constituent)) {
      
      constituent <- unique(as.character(data$constituent))
    }
    
    constituent <- paste( unique(constituent),collapse = ", ")
    
    if (inherits(data, "sf")) {
      data <- sf::st_drop_geometry(data)
    }
    
    fit_res_df <- fit_res_df_st(mod = mod,data = data,idat = idat,...)
    
    ggplot(fit_res_df,aes(x = make_datetime(year_mid))) +
      geom_point(aes(y = y),color = "black",alpha = 0.6) +
      geom_line( aes(y = fit.mean,color = model)) +
      geom_ribbon(aes(ymin = Q1.fit,ymax = Q3.fit,fill = model),alpha = 0.2) +
      facet_grid(st_id ~ model,scales = "free") +
      labs(x = "Year",y = "Value",title = sprintf( "Observed and predicted values [%s]",constituent)) +
      scale_x_datetime(date_breaks = "4 years",date_minor_breaks = "2 years",date_labels = "%Y") +
      theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5 ))
  }


mod_obs_fit_st <- function(mod, data, idat, ...) {
  
  fit_res_df <- fit_res_df_st(mod, data, idat, ...)
  
  ## Remove rows that cannot be plotted
  fit_res_df <- fit_res_df |>
    filter(!is.na(y), !is.na(fit.mean), is.finite(y), is.finite(fit.mean))
  
  ## Use the same limits for observed and fitted values
  plot_limits <- range(c(fit_res_df$y, fit_res_df$fit.mean), na.rm = TRUE)
  
  ## Base plot
  p <- fit_res_df |>
    ggplot(aes(x = y, y = fit.mean, color = st_id)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                alpha = 0.6, col = 2) +
    facet_wrap(~ model, ncol = 3) +
    scale_x_continuous(limits = plot_limits, n.breaks = 4) +
    scale_y_continuous(limits = plot_limits, n.breaks = 4) +
    labs(
      x = "Observed",
      y = "Fitted",
      color = "Station",
      title = sprintf("Observed vs Fitted Values [%s]", const)
    ) +
    guides(color = guide_legend(ncol = 4, byrow = TRUE)) +
    theme(
      plot.title = element_text(size = 14, margin = margin(b = 10)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
  
  ## Extract legend
  legend <- cowplot::get_legend(
    p +
      theme(legend.position = "bottom")
  )
  
  ## Main plot without legend
  p_main <- p +
    theme(legend.position = "none")
  
  ## Convert legend grob to a plot
  p_legend <- cowplot::ggdraw() +
    cowplot::draw_grob(legend)
  
  ## Return separately
  list(
    plot = p_main,
    legend = p_legend
  )
}

mod_obs_res_st <- function(
    mod,
    data,
    idat,
    constituent = NULL,
    ...
) {
  
  ## Get constituent name
  if (is.null(constituent)) {
    
    if ("constituent" %in% names(data)) {
      constituent <- unique(as.character(data$constituent))
    } else {
      constituent <- ""
    }
  }
  
  constituent <- paste(unique(constituent), collapse = ", ")
  
  ## Geometry is not required for this diagnostic plot
  if (inherits(data, "sf")) {
    data <- sf::st_drop_geometry(data)
  }
  
  ## Extract fitted values for every model
  fit_res_df <- fit_res_df_st(
    mod = mod,
    data = data,
    idat = idat,
    ...
  )
  
  ## Ensure the observed-value column is called y
  if (!"y" %in% names(fit_res_df)) {
    
    if ("value" %in% names(fit_res_df)) {
      fit_res_df <- fit_res_df |>
        dplyr::rename(y = value)
    } else {
      stop(
        "`fit_res_df_st()` did not return an observed-value ",
        "column called `y` or `value`."
      )
    }
  }
  
  ## Calculate residuals
  fit_res_df <- fit_res_df |>
    dplyr::mutate(residual = y - fit.mean) |>
    dplyr::filter(
      !is.na(y),
      !is.na(fit.mean),
      !is.na(residual),
      is.finite(y),
      is.finite(fit.mean),
      is.finite(residual)
    )
  
  ## Base residual plot
  p <- fit_res_df |>
    ggplot(aes(x = y, y = residual, color = st_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_point(alpha = 0.6) +
    facet_wrap(~ model, ncol = 3, scales = "free") +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    labs(
      x = "Observed value",
      y = "Residual",
      color = "Station",
      title = sprintf(
        "Observed values versus residuals [%s]",
        constituent
      )
    ) +
    guides(color = guide_legend(ncol = 4, byrow = TRUE)) +
    theme(
      plot.title = element_text(size = 14, margin = margin(b = 10)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    )
  
  ## Extract station legend
  legend <- cowplot::get_legend(
    p +
      theme(
        legend.position = "bottom"
      )
  )
  
  ## Main plot without legend
  p_main <- p +
    theme(
      legend.position = "none"
    )
  
  ## Convert legend grob to a drawable plot
  p_legend <- cowplot::ggdraw() +
    cowplot::draw_grob(legend)
  
  ## Return separately
  list(
    plot = p_main,
    legend = p_legend
  )
}

##==============================================================================
plot_pred_maps <- function(data, fill_var = "mean") {
  lname <- sprintf("Predicted %s [%s]", fill_var, const)
  
  data <- data |>
    mutate(year_mid = as.integer(year_mid))
  
  ggplot() +
    geom_tile(
      data = data,
      aes(x = lon, y = lat, fill = .data[[fill_var]])
    ) +
    geom_sf(data = barrier_km) +
    scale_fill_viridis(name = lname) +
    facet_wrap(~year_mid) +
    theme_void() +
    coord_sf(datum = kmproj)
}


small_legend <- function(
                         direction = "horizontal",
                         barwidth  = 6,
                         barheight = 0.3,
                         title_size = 8,
                         text_size  = 6
                         ) {
    list(
        guides(
            fill = guide_colorbar(
                direction      = direction,
                barwidth       = unit(barwidth, "cm"),
                barheight      = unit(barheight, "cm"),
                title.position = "top",
                title.hjust    = 0.5
            )
        ),
        theme(
            legend.position   = "bottom",
            legend.title      = element_text(size = title_size),
            legend.text       = element_text(size = text_size),
            legend.key.width  = unit(0.3, "cm"),
            legend.key.height = unit(0.3, "cm"),
            legend.box.margin = margin(0, 5, 0, 5)
        )
    )
}



##------------------------------------------------------------------------------
##' @title Calculate Lowest and Highest Astronomical Tide
##'
##' @description Calculates Lowest Astronomical Tide (LAT) and Highest
##'     Astronomical Tide (HAT) from a fitted TideHarmonics model. The function
##'     first predicts the tide hourly over the full astronomical period and
##'     then refines the estimated minimum and maximum using one-minute
##'     predictions around the preliminary extrema.
##'
##' @param fit A fitted tidal model of class \code{tide}, such as an object
##'     returned by \code{fit_th()}.
##' @param from Start of the prediction period. Defaults to
##'     \code{"2000-01-01 00:00:00"}.
##' @param to End of the prediction period. Defaults to
##'     \code{"2018-12-31 23:00:00"}.
##' @param coarse_by Prediction interval, in hours, used for the initial search.
##'     Defaults to 1 hour.
##' @param refine_hours Number of hours on either side of the preliminary
##'     minimum and maximum over which to perform the refined search.
##'     Defaults to 2 hours.
##' @param fine_by Prediction interval, in hours, used for the refined search.
##'     Defaults to 1/60, corresponding to one minute.
##'
##' @return A tibble containing LAT, HAT, their predicted occurrence times,
##'     and the prediction period used.
##'
##' @export
calculate_lat_hat <- function(
    fit,
    from = NULL,
    to = NULL,
    coarse_by = 1,
    refine_hours = 2,
    fine_by = 1/60
) {
  
  ## -------------------------------------------------------------------------
  ## Check fitted object
  ## -------------------------------------------------------------------------
  
  if(!inherits(fit, "tide")) {
    stop(
      "`fit` must be an object of class 'tide', ",
      "such as an object returned by fit_th()."
    )
  }
  
  
  ## -------------------------------------------------------------------------
  ## Check prediction settings
  ## -------------------------------------------------------------------------
  
  if(!is.numeric(coarse_by) ||
     length(coarse_by) != 1L ||
     is.na(coarse_by) ||
     coarse_by <= 0) {
    
    stop(
      "`coarse_by` must be one positive numeric value in hours."
    )
  }
  
  if(!is.numeric(fine_by) ||
     length(fine_by) != 1L ||
     is.na(fine_by) ||
     fine_by <= 0) {
    
    stop(
      "`fine_by` must be one positive numeric value in hours."
    )
  }
  
  if(!is.numeric(refine_hours) ||
     length(refine_hours) != 1L ||
     is.na(refine_hours) ||
     refine_hours <= 0) {
    
    stop(
      "`refine_hours` must be one positive numeric value."
    )
  }
  
  
  ## -------------------------------------------------------------------------
  ## Determine prediction period
  ##
  ## Previously this was fixed at 2000-2018 for every fitted model.
  ## This can fail for newer fitted models because predict.tide() may not
  ## produce valid predictions over a period far outside the model origin.
  ##
  ## Instead, use the fitted model's own origin and search approximately one
  ## complete 19-year lunar nodal cycle.
  ## -------------------------------------------------------------------------
  
  if(is.null(from)) {
    
    if(is.null(fit$orig)) {
      stop(
        "The fitted tide object does not contain `orig`, so a ",
        "prediction start date could not be determined."
      )
    }
    
    fit_origin <- as.POSIXct(
      fit$orig,
      tz = "UTC"
    )
    
    fit_origin <- fit_origin[
      !is.na(fit_origin)
    ]
    
    if(length(fit_origin) == 0L) {
      stop(
        "Could not determine a valid prediction origin from `fit$orig`."
      )
    }
    
    ## Start from 1 January of the model-origin year
    fit_year <- lubridate::year(
      min(
        fit_origin,
        na.rm = TRUE
      )
    )
    
    from <- lubridate::make_datetime(
      year = fit_year,
      month = 1,
      day = 1,
      hour = 0,
      min = 0,
      sec = 0,
      tz = "UTC"
    )
    
  } else {
    
    from <- as.POSIXct(
      from,
      tz = "UTC"
    )
  }
  
  
  if(is.null(to)) {
    
    ## Search approximately one complete 19-year lunar nodal cycle
    to <- from +
      lubridate::years(19) -
      lubridate::hours(1)
    
  } else {
    
    to <- as.POSIXct(
      to,
      tz = "UTC"
    )
  }
  
  
  ## -------------------------------------------------------------------------
  ## Check prediction limits
  ## -------------------------------------------------------------------------
  
  if(length(from) != 1L ||
     length(to) != 1L ||
     is.na(from) ||
     is.na(to)) {
    
    stop(
      "`from` and `to` must be valid single date-time values."
    )
  }
  
  if(from >= to) {
    stop(
      "`from` must be earlier than `to`."
    )
  }
  
  
  ## -------------------------------------------------------------------------
  ## Coarse search
  ##
  ## predict.tide() uses hours for the `by` argument.
  ## An hourly interval therefore corresponds to by = 1.
  ## -------------------------------------------------------------------------
  
  pred_coarse <- predict(
    fit,
    from = from,
    to = to,
    by = coarse_by
  )
  
  pred_coarse <- as.numeric(
    pred_coarse
  )
  
  
  ## -------------------------------------------------------------------------
  ## Check coarse predictions
  ## -------------------------------------------------------------------------
  
  if(length(pred_coarse) == 0L ||
     all(is.na(pred_coarse))) {
    
    stop(
      paste0(
        "No valid tidal predictions were produced in the coarse search. ",
        "Prediction period: ",
        format(from),
        " to ",
        format(to),
        "."
      )
    )
  }
  
  
  ## Remove non-finite predictions if any are present
  valid_coarse <- is.finite(
    pred_coarse
  )
  
  if(!any(valid_coarse)) {
    stop(
      paste0(
        "No finite tidal predictions were produced in the coarse search. ",
        "Prediction period: ",
        format(from),
        " to ",
        format(to),
        "."
      )
    )
  }
  
  
  ## -------------------------------------------------------------------------
  ## Reconstruct the times associated with the predictions
  ## -------------------------------------------------------------------------
  
  time_coarse <- seq(
    from = from,
    by = coarse_by * 60 * 60,
    length.out = length(pred_coarse)
  )
  
  
  ## Find HAT/LAT using finite predictions only
  pred_coarse_valid <- pred_coarse[
    valid_coarse
  ]
  
  time_coarse_valid <- time_coarse[
    valid_coarse
  ]
  
  idx_hat_coarse <- which.max(
    pred_coarse_valid
  )
  
  idx_lat_coarse <- which.min(
    pred_coarse_valid
  )
  
  time_hat_coarse <- time_coarse_valid[
    idx_hat_coarse
  ]
  
  time_lat_coarse <- time_coarse_valid[
    idx_lat_coarse
  ]
  
  
  ## -------------------------------------------------------------------------
  ## Fine search around the preliminary HAT
  ## -------------------------------------------------------------------------
  
  time_hat_from <-
    time_hat_coarse -
    refine_hours * 60 * 60
  
  time_hat_to <-
    time_hat_coarse +
    refine_hours * 60 * 60
  
  
  pred_hat_fine <- predict(
    fit,
    from = time_hat_from,
    to = time_hat_to,
    by = fine_by
  )
  
  pred_hat_fine <- as.numeric(
    pred_hat_fine
  )
  
  
  if(length(pred_hat_fine) == 0L ||
     all(is.na(pred_hat_fine))) {
    
    stop(
      "No valid tidal predictions were produced around the estimated HAT."
    )
  }
  
  
  valid_hat <- is.finite(
    pred_hat_fine
  )
  
  if(!any(valid_hat)) {
    stop(
      "No finite tidal predictions were produced around the estimated HAT."
    )
  }
  
  
  time_hat_fine <- seq(
    from = time_hat_from,
    by = fine_by * 60 * 60,
    length.out = length(pred_hat_fine)
  )
  
  pred_hat_valid <- pred_hat_fine[
    valid_hat
  ]
  
  time_hat_valid <- time_hat_fine[
    valid_hat
  ]
  
  idx_hat <- which.max(
    pred_hat_valid
  )
  
  
  ## -------------------------------------------------------------------------
  ## Fine search around the preliminary LAT
  ## -------------------------------------------------------------------------
  
  time_lat_from <-
    time_lat_coarse -
    refine_hours * 60 * 60
  
  time_lat_to <-
    time_lat_coarse +
    refine_hours * 60 * 60
  
  
  pred_lat_fine <- predict(
    fit,
    from = time_lat_from,
    to = time_lat_to,
    by = fine_by
  )
  
  pred_lat_fine <- as.numeric(
    pred_lat_fine
  )
  
  
  if(length(pred_lat_fine) == 0L ||
     all(is.na(pred_lat_fine))) {
    
    stop(
      "No valid tidal predictions were produced around the estimated LAT."
    )
  }
  
  
  valid_lat <- is.finite(
    pred_lat_fine
  )
  
  if(!any(valid_lat)) {
    stop(
      "No finite tidal predictions were produced around the estimated LAT."
    )
  }
  
  
  time_lat_fine <- seq(
    from = time_lat_from,
    by = fine_by * 60 * 60,
    length.out = length(pred_lat_fine)
  )
  
  pred_lat_valid <- pred_lat_fine[
    valid_lat
  ]
  
  time_lat_valid <- time_lat_fine[
    valid_lat
  ]
  
  idx_lat <- which.min(
    pred_lat_valid
  )
  
  
  ## -------------------------------------------------------------------------
  ## Return results
  ## -------------------------------------------------------------------------
  
  out <- tibble(
    constituent = c(
      "LAT",
      "HAT"
    ),
    
    value = c(
      pred_lat_valid[idx_lat],
      pred_hat_valid[idx_hat]
    ),
    
    time = c(
      time_lat_valid[idx_lat],
      time_hat_valid[idx_hat]
    ),
    
    prediction_from = from,
    prediction_to = to
  )
  
  
  return(out)
}


show_temporal_results <- function(
    constituent_name,
    res,
    db,
    min_obs = MIN_TEMPORAL_OBS
) {
  
  ## Models for this constituent
  models <- res[[constituent_name]]
  
  ## Number of usable temporal observations
  n_obs <- sum(
    !is.na(
      db[[constituent_name]]$value
    )
  )
  
  
  ## -------------------------------------------------------------------------
  ## No temporal models available
  ## -------------------------------------------------------------------------
  
  if (
    is.null(models) ||
    length(models) == 0
  ) {
    
    cat(
      sprintf(
        paste0(
          "\n\n**Temporal modelling was not performed for %s.** ",
          "Only %i non-missing temporal estimates were available; ",
          "at least %i are required for temporal modelling.\n\n"
        ),
        constituent_name,
        n_obs,
        min_obs
      )
    )
    
    return(
      invisible(NULL)
    )
  }
  
  
  ## -------------------------------------------------------------------------
  ## Model summary table
  ## -------------------------------------------------------------------------
  
  mod_table <- mod_summ(
    models,
    db[[constituent_name]]
  )
  
  print(
    out_table(
      mod_table
    )
  )
  
  
  ## -------------------------------------------------------------------------
  ## Posterior marginals
  ## -------------------------------------------------------------------------
  
  p_marg <- mod_marg(
    models,
    constituent = constituent_name
  )
  
  print(
    p_marg
  )
  
  
  ## -------------------------------------------------------------------------
  ## Observed and fitted values over time
  ## -------------------------------------------------------------------------
  
  p_pred <- mod_pred(
    models,
    db[[constituent_name]],
    constituent = constituent_name
  )
  
  print(
    p_pred
  )
  
  
  ## -------------------------------------------------------------------------
  ## Observed versus fitted
  ## -------------------------------------------------------------------------
  
  p_obs_fit <- mod_obs_fit(
    models,
    db[[constituent_name]],
    constituent = constituent_name
  )
  
  print(
    p_obs_fit
  )
  
  
  ## -------------------------------------------------------------------------
  ## Observed versus residuals
  ## -------------------------------------------------------------------------
  
  p_obs_res <- mod_obs_res(
    models,
    db[[constituent_name]],
    constituent = constituent_name
  )
  
  print(
    p_obs_res
  )
  
  
  invisible(NULL)
}

##------------------------------------------------------------------------------
## Create spatial prediction plot only when it does not already exist

save_pred_plot <- function(const, fill_var) {
  
  fl <- here::here("figs", sprintf("pred-%s-%s.png", fill_var, const))
  
  if (!file.exists(fl)) {
    
    cat(sprintf("Creating prediction plot: %s\n", fl))
    
    plot_data <- df_pred[[const]]
  
    
    p <- plot_pred_maps(plot_data, fill_var = fill_var) +
      small_legend()
    
    ggsave(fl, plot = p, width = 12, height = 12, dpi = 100)
    
  } else {
    
    cat(sprintf("Prediction plot already exists, using: %s\n", fl))
  }
  
  fl
}