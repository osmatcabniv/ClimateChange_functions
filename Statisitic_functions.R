# Load all necessary libraries
library(data.table)
library(SPEI)
library(trend)
library(gamlss)
library(gamlss.dist)
library(extRemes) # For fevd() function
library(e1071)    # For skewness() function

# --- Reusable Functions for Seasonal Analysis ---

#' @title Define Precipitation and Temperature Seasons
#' @description Creates a master definition table for a 365-day calendar based on precipitation (wet/dry) and temperature (winter, spring, summer, fall).
#' @param dt A melted data.table with historical data, containing 'Date', 'variable', and 'value'.
#' @param smoothing_window The number of days for the temperature smoothing. Default is 21.
#' @return A data.table with columns: 'day_365', 'pr_season', and 'tmp_season'.
define_seasons <- function(dt, smoothing_window = 21) {
  setDT(dt)
  
  # --- 1. Precipitation Season (Wet/Dry) Definition (Optimized) ---
  # Create a filtered data.table for daily precipitation once to avoid redundant filtering.
  pr_data <- dt[variable == 'pr' & !(month(Date) == 2 & mday(Date) == 29) & !is.na(value)]
  
  # Calculate the overall mean from the daily values.
  overall_mean_pr <- mean(pr_data$value)
  
  # Calculate daily climatology and cumulative anomalies in a single chain.
  daily_pr_clim <- pr_data[, .(mean_daily_pr = mean(value)), by = .(day_365 = rowid(year(Date)))
  ][, anomaly := mean_daily_pr - ..overall_mean_pr
  ][, cum_anomaly := cumsum(anomaly)]
  setkey(daily_pr_clim, day_365)
  
  wet_start <- daily_pr_clim[which.min(cum_anomaly), day_365] + 1
  wet_end <- daily_pr_clim[which.max(cum_anomaly), day_365]
  
  pr_season_def <- data.table(day_365 = 1:365)
  if (wet_start > wet_end) {
    pr_season_def[, pr_season := ifelse(day_365 >= wet_start | day_365 <= wet_end, "wet", "dry")]
  } else {
    pr_season_def[, pr_season := ifelse(day_365 >= wet_start & day_365 <= wet_end, "wet", "dry")]
  }
  
  # --- 2. Temperature Season Definition (already optimized) ---
  daily_clim <- dt[variable %in% c("tasmin", "tasmax")
  ][, dcast(.SD, Date ~ variable, value.var = "value")
  ][, tavg := (tasmax + tasmin) / 2
  ][!is.na(tavg) & !(month(Date) == 2 & mday(Date) == 29),
    .(tavg, day_365 = rowid(year(Date)))
  ][, .(mean_tavg = mean(tavg)), by = day_365]
  
  pad <- floor(smoothing_window / 2)
  extended_tavg <- c(tail(daily_clim$mean_tavg, pad), daily_clim$mean_tavg, head(daily_clim$mean_tavg, pad))
  daily_clim$smoothed_tavg <- frollmean(extended_tavg, n = smoothing_window, align = "center")[(pad + 1):(pad + 365)]
  
  quantiles <- quantile(daily_clim$mean_tavg, probs = c(0.25, 0.75), na.rm = TRUE)
  daily_clim[smoothed_tavg <= quantiles[1], season := "winter"]
  daily_clim[smoothed_tavg >= quantiles[2], season := "summer"]
  
  find_main_season_block <- function(seasons_vec, target_season) {
    rle_obj <- rle(c(seasons_vec, seasons_vec))
    runs <- which(rle_obj$values == target_season)
    if (length(runs) == 0) return(list(start = -1, end = -1))
    main_run <- runs[which.max(rle_obj$lengths[runs])]
    run_end_730 <- cumsum(rle_obj$lengths)[main_run]
    run_start_730 <- run_end_730 - rle_obj$lengths[main_run] + 1
    return(list(start = (run_start_730 - 1) %% 365 + 1, end = (run_end_730 - 1) %% 365 + 1))
  }
  w <- find_main_season_block(daily_clim$season, "winter")
  s <- find_main_season_block(daily_clim$season, "summer")
  
  daily_clim[, final_season := "transition"]
  if (w$start > w$end) { daily_clim[day_365 >= w$start | day_365 <= w$end, final_season := "winter"] } else { daily_clim[day_365 >= w$start & day_365 <= w$end, final_season := "winter"] }
  if (s$start > s$end) { daily_clim[day_365 >= s$start | day_365 <= s$end, final_season := "summer"] } else { daily_clim[day_365 >= s$start & day_365 <= s$end, final_season := "summer"] }
  
  if (s$start < w$start) { # Northern Hemisphere
    daily_clim[final_season == "transition" & day_365 > w$end & day_365 < s$start, final_season := "spring"]
    daily_clim[final_season == "transition", final_season := "fall"]
  } else { # Southern Hemisphere
    daily_clim[final_season == "transition" & day_365 > s$end & day_365 < w$start, final_season := "fall"]
    daily_clim[final_season == "transition", final_season := "spring"]
  }
  
  tmp_season_def <- daily_clim[, .(day_365, tmp_season = final_season)]
  
  # --- 3. Combine into Master Definition Table ---
  master_seasons <- pr_season_def[tmp_season_def, on = "day_365"]
  setnames(master_seasons, c("DAY", "PR_SEASON", "TMP_SEASON")) # Use standard names
  return(master_seasons)
}

#' Gets return levels from the best GEV/Gumbel model using a true worst-case scenario.
#'
#' This function automatically tests for trends, selects the best model, and
#' computes return levels. For non-stationary models, it calculates levels for
#' both the start and end of the period, compares them, and returns the set of
#' values corresponding to the highest risk.
#'
#' @param x A numeric vector of annual maxima.
#' @param return_periods A numeric vector of return periods.
#' @param alpha The significance level for statistical tests.
#' @return A data.frame containing the calculated true worst-case return levels.
get_true_worst_case_return_levels <- function(x, return_periods = c(5, 10, 25, 50, 100, 200, 500), alpha = 0.05) {
  
  mk_test_result <- trend::mk.test(x)
  
  # --- Non-Stationary Path ---
  if (mk_test_result$p.value < alpha) {
    # --- 1. Non-Stationary Model Selection ---
    model_structures <- list(
      `Location-NS` = list(location.fun = ~ time, type = "GEV"),
      `Scale-NS` = list(scale.fun = ~ time, type = "GEV"), 
      `Location_Scale-NS` = list(location.fun = ~ time, scale.fun = ~ time, type = "GEV")
    )
    time <- 1:length(x)
    fit_results <- lapply(model_structures, function(args) {
      fit <- try(do.call("fevd", c(list(x = x, time = time), args)), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      return(list(fit = fit, AIC = 2 * length(fit$results$par) + 2 * fit$results$value))
    })
    successful_fits <- fit_results[!sapply(fit_results, is.null)]
    if (length(successful_fits) == 0) stop("All non-stationary models failed to converge.")
    aics <- sapply(successful_fits, `[[`, "AIC")
    final_model_nonstat <- successful_fits[[which.min(aics)]]$fit
    
    # --- 2. Calculate Non-Stationary Worst-Case Return Levels ---
    time_points <- c(1, length(x))
    rl_matrix_nonstat <- return.level(final_model_nonstat, return.period = return_periods, t_covars = data.frame(time = time_points))
    sen_results <- trend::sens.slope(x)
    
    if (rl_matrix_nonstat[1, ncol(rl_matrix_nonstat)] >= rl_matrix_nonstat[2, ncol(rl_matrix_nonstat)]) {
      worst_case_row <- rl_matrix_nonstat[1, ]
      time_point_used <- 1
      scenario <- "Start of Period (Highest Risk)"
    } else {
      worst_case_row <- rl_matrix_nonstat[length(x), ]
      time_point_used <- length(x)
      scenario <- "End of Period (Highest Risk)"
    }
    
    rl_df_nonstat <- as.data.frame(t(worst_case_row))
    colnames(rl_df_nonstat) <- names(worst_case_row)
    rl_df_nonstat$`MK p-value` <- mk_test_result$p.value
    rl_df_nonstat$`Sen p-value` <- sen_results$p.value
    rl_df_nonstat$`Sen slope` <- sen_results$estimates
    rl_df_nonstat$Worst_Case_Scenario <- scenario
    rl_df_nonstat$Time_Step_Used <- time_point_used
    
    # --- 3. Calculate Stationary Return Levels for Comparison ---
    stat_fit <- fevd(x, type = "GEV")
    shape_se_stat <- sqrt(diag(solve(stat_fit$results$hessian)))["shape"]
    if (!is.nan(shape_se_stat) && abs(stat_fit$results$par["shape"] / shape_se_stat) < qnorm(1 - alpha / 2)) {
      final_model_stat <- fevd(x, type = "Gumbel")
    } else {
      final_model_stat <- stat_fit
    }
    rl_matrix_stat <- return.level(final_model_stat, return.period = return_periods)
    rl_df_stat <- as.data.frame(t(rl_matrix_stat))
    rl_df_stat <- as.data.frame(t(rl_df_stat$x[1,]))
    colnames(rl_df_stat) <- names(worst_case_row)
    
    # --- 4. Combine and return ---
    return(rbindlist(list(rl_df_nonstat, rl_df_stat), use.names = TRUE, fill = TRUE))
    
  } else {
    # --- Stationary Path (No Trend Detected) ---
    final_model <- fevd(x, type = "GEV")
    
    rl_matrix <- return.level(final_model, return.period = return_periods)
    rl_df <- as.data.frame(t(rl_matrix))
    rl_df <- as.data.frame(t(rl_df$x[1,]))
    colnames(rl_df) <- names(rl_matrix)
    rl_df$`MK p-value` <- mk_test_result$p.value
    rl_df$Worst_Case_Scenario <- "Stationary"
    return(rl_df)
  }
}

#' Calculate SPI using the highly optimized SPEI package.
#'
#' @param dt A data.table with 'year', 'month', and 'TotalPrecip' columns.
#' @param scales A numeric vector of scales (e.g., c(3, 6, 12)).
#' @return A data.table with new SPI columns (e.g., SPI_3, SPI_6).
#' 
calculate_spi_package <- function(dt, scales = c(3,6,12)) {
  dt_spi <- copy(dt)
  start_year <- dt_spi[1, year]
  start_month <- dt_spi[1, month]
  precip_ts <- ts(dt_spi$TotalPrecip, start = c(start_year, start_month), frequency = 12)
  for (s in scales) {
    spi_col <- paste0("SPI_", s)
    spi_result <- spi(precip_ts, scale = s, na.rm = TRUE, verbose = F)
    dt_spi[, (spi_col) := as.numeric(spi_result$fitted)]
  }
  return(dt_spi)
}

#' Optimized calculation and comparison of stationary and nonstationary SPI.
#'
#' @param dt A data.table with 'year', 'month', and 'TotalPrecip' columns.
#' @param scales A numeric vector of scales (e.g., c(3, 6, 12)).
#' @return A list containing 'results' and 'summary'.

calculate_spi_analysis_optimized <- function(dt, scales = c(3,6,12), alpha = 0.05) {
  
  # --- STEP 1: EFFICIENT PRE-COMPUTATION ---
  
  # Calculate all stationary SPIs and moving sums in one go.
  #message("Step 1: Pre-computing all stationary SPIs and moving sums...")
  dt_final <- calculate_spi_package(dt, scales)
  for (s in scales) {
    sum_col <- paste0("P_", s)
    nspi_col <- paste0("NSPI_", s)
    dt_final[, (sum_col) := frollsum(TotalPrecip, n = s, align = "right")]
    dt_final[, (nspi_col) := .SD[[paste0("SPI_", s)]]] # Create NSPI cols
  }
  
  # --- STEP 2: BATCH TREND ANALYSIS USING data.table's 'by' ---
  
  #message("Step 2: Performing all trend tests in a single grouped operation...")
  
  # Melt the data from wide to long format for efficient grouping
  measure_cols <- paste0("P_", scales)
  dt_long <- melt(dt_final, 
                  id.vars = "month", 
                  measure.vars = measure_cols, 
                  variable.name = "scale_col", 
                  value.name = "prec_sum",
                  variable.factor = FALSE)
  
  # Run trend tests for all months and scales at once. This is much faster than an R loop.
  trend_summary <- dt_long[!is.na(prec_sum), {
    mk_test <- trend::mk.test(prec_sum)
    sen_test <- trend::sens.slope(prec_sum)
    .(
      `MK p-value` = mk_test$p.value,
      `Sen p-value` = sen_test$p.value,
      `Sen slope` = sen_test$estimates
    )
  }, by = .(scale_col, month)]
  
  # --- STEP 3: TARGETED NONSTATIONARY FITTING ---
  
  #message("Step 3: Fitting nonstationary models only where significant trends were found...")
  
  # Identify only the month/scale combos that need the expensive GAMLSS fit
  trends_to_fix <- trend_summary[`MK p-value` < alpha]
  
  if (nrow(trends_to_fix) > 0) {
    # Loop through the much smaller list of detected trends
    for (i in 1:nrow(trends_to_fix)) {
      current_scale_col <- trends_to_fix$scale_col[i]
      current_month <- trends_to_fix$month[i]
      
      s <- as.integer(gsub("P_", "", current_scale_col))
      nspi_col <- paste0("NSPI_", s)
      
      #message(paste0("  -> Fitting NSPI for scale ", s, ", month ", current_month, "..."))
      
      # Extract data for fitting
      monthly_sums <- dt_final[month == current_month, get(current_scale_col)]
      fit_data <- data.frame(prec = monthly_sums, time = 1:length(monthly_sums))
      fit_data <- na.omit(fit_data)
      
      # Fit the nonstationary model
      nonstat_fit <- try(
        gamlss(prec ~ time, family = ZAGA(), data = fit_data, trace = FALSE),
        silent = TRUE
      )
      
      if (!inherits(nonstat_fit, "try-error")) {
        cdf_vals <- pZAGA(fit_data$prec,
                          mu = predict(nonstat_fit, what = "mu", type = "response"),
                          sigma = predict(nonstat_fit, what = "sigma", type = "response"),
                          nu = predict(nonstat_fit, what = "nu", type = "response"))
        
        nonstat_spi <- qnorm(cdf_vals)
        
        # Update the NSPI column in the main results table
        dt_final[month == current_month & !is.na(get(current_scale_col)), (nspi_col) := nonstat_spi]
      }
    }
  } else {
    #message("  -> No significant trends found at any scale.")
  }
  
  # --- STEP 4: FINAL CLEANUP AND FORMATTING ---
  
  #message("Step 4: Finalizing output...")
  
  # Clean up intermediate P_ columns
  dt_final[, (measure_cols) := NULL]
  
  # Format the summary list
  summary_list <- split(trend_summary, by = "scale_col")
  summary_list <- lapply(summary_list, function(df) {
    df[, scale_col := NULL]
    setorder(df, month)
    df[, Month := month.abb[month]][, month := NULL]
    setcolorder(df, "Month")
    return(df)
  })
  names(summary_list) <- gsub("P_", "scale_", names(summary_list))
  
  # Remove NSPI columns where no trends were found
  scales_with_trends <- unique(as.integer(gsub("P_", "", trends_to_fix$scale_col)))
  scales_without_trends <- setdiff(scales, scales_with_trends)
  
  if (length(scales_without_trends) > 0) {
    for (s in scales_without_trends) {
      nspi_col_to_remove <- paste0("NSPI_", s)
      #message(paste0("No trends for scale ", s, ". Removing ", nspi_col_to_remove, "."))
      dt_final[, (nspi_col_to_remove) := NULL]
    }
  }
  
  #message("Analysis complete.")
  return(list(results = dt_final, summary = summary_list))
}


#' @title Calculate All Statistics Using a Pre-defined Season Map (Optimized)
#' @description Calculates seasonal, monthly, and yearly statistics for climate data. This version is optimized for performance.
#' @param dt A melted data.table with climate data.
#' @param season_map A data.table with season definitions, from the `define_seasons` function.
#' @param calculate_daily_avg A logical flag. If TRUE, computes the 365-day climatology.
#' @param report_yearly_summary A logical flag. If TRUE, computes a detailed yearly summary time series.
#' @return A list of data.tables with all the calculated statistics.
calculate_statistics <- function(dt, season_map, calculate_daily_avg = FALSE, report_yearly_summary = FALSE, scales = c(3,6,12), return_periods = c(5, 10, 25, 50, 100, 200, 500), alpha = 0.05) {
  
  setDT(dt)
  setnames(season_map, c("day_365", "pr_season", "tmp_season"))
  
  # --- 1. Prepare a unified data table ---
  pr_dt <- dt[variable == 'pr']
  tavg_dt <- dcast(dt[variable %in% c('tasmin', 'tasmax')], 
                   Date + source + year ~ variable, 
                   value.var = "value"
  )[, .(Date, variable = "tavg", source, value = (tasmin + tasmax) / 2, year)
  ][!is.na(value)]
  
  # Rbind once at the beginning
  analysis_dt <- rbind(pr_dt, tavg_dt)
  
  # --- 2. Perform a single join with the season map ---
  analysis_dt[, month := month(Date)]
  dt_seasonal <- season_map[analysis_dt[!(month == 2 & mday(Date) == 29)][, day_365 := rowid(year), by = .(variable)], on = "day_365", nomatch = 0]
  
  # --- 3. Define a single, powerful aggregation function ---
  calculate_summary <- function(data, time_vars) {
    agg_yearly <- data[, .(agg_value = if(variable[1] == 'pr') sum(value) else mean(value)), 
                       by = .(variable, year, get(time_vars))]
    setnames(agg_yearly, "get", time_vars)
    
    summary <- agg_yearly[, {
      quants <- quantile(agg_value, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      .(mean = mean(agg_value), 
        std = sd(agg_value), 
        skewness = skewness(agg_value),
        min = min(agg_value),
        Q1 = quants[1],
        Q2 = quants[2],
        Q3 = quants[3],
        max = max(agg_value))
    }, by = .(variable, get(time_vars))]
    setnames(summary, "get", time_vars)
    return(summary)
  }
  
  # --- 4. Calculate all statistics using the unified function ---
  by_pr_season  <- calculate_summary(dt_seasonal, "pr_season")
  by_tmp_season <- calculate_summary(dt_seasonal, "tmp_season")
  monthly       <- calculate_summary(dt_seasonal, "month")
  
  # --- 5. Calculate Yearly Time Series directly ---
  yearly_timeseries <- dt_seasonal[, .(value = if(variable[1] == 'pr') sum(value) else mean(value)), by = .(variable, year)]
  setkey(yearly_timeseries, variable, year)
  
  # --- 6. Create annual summary for monthly table ---
  annual_summary <- yearly_timeseries[, {
    quants <- quantile(value, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    .(mean = mean(value), 
      std = sd(value), 
      skewness = skewness(value),
      min = min(value),
      Q1 = quants[1],
      Q2 = quants[2],
      Q3 = quants[3],
      max = max(value))
  }, by = variable][, month := 13]
  
  final_monthly <- rbind(monthly, annual_summary, fill = TRUE)
  setkey(final_monthly, variable, month)
  
  # --- 7. Prepare final result list ---
  results_list <- list(
    by_pr_season = by_pr_season,
    by_tmp_season = by_tmp_season,
    monthly = final_monthly,
    yearly = yearly_timeseries
  )
  
  # --- 8. Conditionally calculate daily average ---
  if (calculate_daily_avg) {
    daily_average <- dt_seasonal[, .(mean_value = mean(value)), by = .(variable, day_365)]
    setkey(daily_average, variable, day_365)
    results_list$daily_average <- daily_average
  }
  
  # --- 9. Calculate GEV return levels for pr and tavg ---
  gev_results_list <- list()
  # For precipitation
  annual_maxima_pr <- dt_seasonal[variable == 'pr', .(max_val = max(value)), by = year]
  if(nrow(annual_maxima_pr) > 0 && length(unique(annual_maxima_pr$max_val)) > 1) { 
    gev_results_list$pr <- get_true_worst_case_return_levels(annual_maxima_pr$max_val, return_periods, alpha)
  }
  # For average temperature
  annual_maxima_tavg <- dt_seasonal[variable == 'tavg', .(max_val = max(value)), by = year]
  if(nrow(annual_maxima_tavg) > 0 && length(unique(annual_maxima_tavg$max_val)) > 1) { 
    gev_results_list$tavg <- get_true_worst_case_return_levels(annual_maxima_tavg$max_val, return_periods, alpha)
  }
  if(length(gev_results_list) > 0){
    results_list$return_levels <- gev_results_list
  }
  
  # --- 10. Calculate SPI for precipitation ---
  monthly_precip_ts <- dt_seasonal[variable == 'pr', .(TotalPrecip = sum(value)), by = .(year, month)]
  setorder(monthly_precip_ts, year, month)
  if(nrow(monthly_precip_ts) > 0) {
    results_list$spi <- calculate_spi_analysis_optimized(monthly_precip_ts, scales, alpha)
  }
  
  # --- 11. Conditionally generate detailed yearly summary time series ---
  if (report_yearly_summary) {
    # Monthly data in wide format
    monthly_agg <- dt_seasonal[, .(value = if(variable[1] == 'pr') sum(value) else mean(value)), by = .(variable, year, month)]
    monthly_wide <- dcast(monthly_agg, variable + year ~ month, value.var = "value")
    setnames(monthly_wide, as.character(1:12), month.abb)
    
    # Yearly total/average
    yearly_total <- copy(yearly_timeseries)
    setnames(yearly_total, "value", "Yearly")
    
    # Seasonal aggregations
    pr_season_wide <- dcast(dt_seasonal[, .(value = if(variable[1] == 'pr') sum(value) else mean(value)), by = .(variable, year, pr_season)],
                            variable + year ~ pr_season, value.var = "value")
    tmp_season_wide <- dcast(dt_seasonal[, .(value = if(variable[1] == 'pr') sum(value) else mean(value)), by = .(variable, year, tmp_season)],
                             variable + year ~ tmp_season, value.var = "value")
    
    # Annual maxima
    annual_maxima <- dt_seasonal[, .(Maxima = max(value)), by = .(variable, year)]
    
    # Join all pieces together
    yearly_summary <- monthly_wide[yearly_total, on = c("variable", "year")
    ][pr_season_wide, on = c("variable", "year")
    ][tmp_season_wide, on = c("variable", "year")
    ][annual_maxima, on = c("variable", "year")]
    
    # Set column order for clarity
    setcolorder(yearly_summary, c("variable", "year", month.abb, "Yearly", "wet", "dry", "winter", "spring", "summer", "fall", "Maxima"))
    
    results_list$yearly_summary <- yearly_summary
    
    # --- 12. Calculate and add trend summary for ALL yearly summary variables ---
    # Melt the summary table to a long format for efficient trend analysis
    cols_to_analyze <- setdiff(names(yearly_summary), c("variable", "year"))
    long_summary <- melt(yearly_summary, 
                         id.vars = c("variable", "year"), 
                         measure.vars = cols_to_analyze, 
                         variable.name = "metric", 
                         value.name = "value")
    
    # Calculate trends for each variable and metric
    trend_summary <- long_summary[!is.na(value), {
      res <- .(MK_p_value = NA_character_, Sen_p_value = NA_character_, Sen_slope = NA_character_)
      # Use tryCatch for maximum robustness in batch processing
      tryCatch({
        if (.N > 3 && length(unique(value)) > 1) {
          mk <- trend::mk.test(value)
          sen <- trend::sens.slope(value)
          res <- .(
            MK_p_value = formatC(mk$p.value, format = "e", digits = 2),
            Sen_p_value = formatC(sen$p.value, format = "e", digits = 2),
            Sen_slope = formatC(sen$estimates, format = "e", digits = 2)
          )
        }
      }, error = function(e) {
        # If an error occurs, 'res' remains as NAs
        warning(paste("Trend analysis failed for variable:", .BY$variable, "and metric:", .BY$metric))
      })
      res
    }, by = .(variable, metric)]
    
    setorder(trend_summary, variable, metric)
    results_list$trend_summary <- trend_summary
    
  }
  
  return(results_list)
}