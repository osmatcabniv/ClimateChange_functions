#=================================================================================#
# SCRIPT TO CALCULATE ETCCDI CLIMATE INDICES EFFICIENTLY (ALL-IN-ONE FUNCTION)
#=================================================================================#
#
# DESCRIPTION:
# This script provides a single, self-contained function `calculate_all_etccdi`
# for calculating the core ETCCDI climate indices. It integrates high-performance
# data preparation using 'data.table' and the final index calculations.
#
# HOW TO USE:
#   1. Save this entire script as an .R file (e.g., "etccdi_calculator.R").
#   2. In your main analysis script, load it using: source("etccdi_calculator.R")
#   3. Call the `calculate_all_etccdi()` function with your data.
#
#=================================================================================#
# SCRIPT SETUP: LOAD REQUIRED PACKAGES
#=================================================================================#
required_packages <- c("data.table", "PCICt", "climdex.pcic", "lubridate")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#=================================================================================#
# MAIN FUNCTION: calculate_all_etccdi
#=================================================================================#
# This function takes daily climate data, prepares it, and calculates indices.
calculate_all_etccdi <- function(
    tmax = NULL, tmin = NULL, prec = NULL,
    tmax.dates = NULL, tmin.dates = NULL, prec.dates = NULL,
    base.range = c(1961, 1990), n = 5,
    northern.hemisphere = TRUE, tavg = NULL, tavg.dates = NULL,
    quantiles = NULL, temp.qtiles = c(0.1, 0.9),
    prec.qtiles = c(0.95, 0.99),
    max.missing.days = c(annual = 15, monthly = 3),
    min.base.data.fraction.present = 0.1,
    calculate_indices = TRUE
) {
  
  #-------------------------------------------------------------------------------#
  # PART 1: DATA PREPARATION
  #-------------------------------------------------------------------------------#
  
  # --- 1. Argument validation ---
  climdex.pcic:::check.basic.argument.validity(tmax, tmin, prec, tmax.dates, tmin.dates, prec.dates, base.range, n, tavg, tavg.dates)
  stopifnot(length(max.missing.days) == 2 && all(c("annual", "monthly") %in% names(max.missing.days)))
  stopifnot(is.numeric(min.base.data.fraction.present) && length(min.base.data.fraction.present) == 1)
  
  # --- 2. Consolidate all data into a single data.table ---
  all_dates_list <- list(tmax.dates, tmin.dates, prec.dates, tavg.dates)
  all_dates_list <- all_dates_list[!sapply(all_dates_list, is.null)]
  if (length(all_dates_list) == 0) stop("At least one of tmax, tmin, prec, or tavg must be provided.")
  
  full_date_range <- range(do.call(c, all_dates_list))
  calendar <- attr(all_dates_list[[1]], "cal")
  
  date_series <- seq(full_date_range[1], full_date_range[2], by = "day")
  clim_data <- data.table(date = date_series)
  
  if (!is.null(tmax)) clim_data[data.table(date = tmax.dates, tmax = tmax), on = "date", tmax := i.tmax]
  if (!is.null(tmin)) clim_data[data.table(date = tmin.dates, tmin = tmin), on = "date", tmin := i.tmin]
  if (!is.null(prec)) clim_data[data.table(date = prec.dates, prec = prec), on = "date", prec := i.prec]
  if (!is.null(tavg)) clim_data[data.table(date = tavg.dates, tavg = tavg), on = "date", tavg := i.tavg]
  
  if (is.null(tavg) && !is.null(tmin) && !is.null(tmax)) {
    clim_data[, tavg := (tmax + tmin) / 2]
  }
  
  # --- 3. Add date components ---
  clim_data[, `:=`(year = year(date), month = month(date), yday = yday(date))]
  
  if(calendar != "365_day") {
    clim_data[month == 2 & mday(date) == 29, yday := NA]
    clim_data[, yday := nafill(yday, type = "locf")]
  }
  
  date_factors <- list(
    annual = factor(clim_data$year),
    monthly = factor(format(clim_data$date, format = "%Y-%m", tz = "GMT"))
  )
  
  # --- 4. Calculate NA masks ---
  var_list <- c("tmax", "tmin", "prec", "tavg")
  present_vars <- var_list[sapply(var_list, function(x) x %in% names(clim_data))]
  
  namasks <- list(monthly = list(), annual = list())
  for (var in present_vars) {
    monthly_na <- clim_data[, .(na_count = sum(is.na(get(var)))), by = .(year, month)]
    monthly_mask <- monthly_na$na_count <= max.missing.days["monthly"]
    
    annual_na <- clim_data[, .(na_count = sum(is.na(get(var)))), by = year]
    annual_mask <- annual_na$na_count <= max.missing.days["annual"]
    
    combined_annual_mask <- as.numeric(tapply(monthly_mask, monthly_na$year, prod))
    final_annual_mask <- pmin(annual_mask, combined_annual_mask)
    
    namasks$monthly[[var]] <- as.numeric(monthly_mask)
    namasks$annual[[var]] <- as.numeric(final_annual_mask)
  }
  
  # --- 5. Calculate Quantiles ---
  quantiles_env <- new.env(parent = emptyenv())
  filled.list <- lapply(clim_data[, ..present_vars], as.numeric)
  
  bs.date.range <- PCICt::as.PCICt(paste(base.range, c("01-01", "12-31"), sep = "-"), cal = calendar)
  
  bs.date.series <- seq(bs.date.range[1], bs.date.range[2], by = "day")
  days_threshold <- 359
  
  if ("tmax" %in% present_vars && length(intersect(tmax.dates, bs.date.series)) > days_threshold) {
    delayedAssign("tmax", climdex.pcic:::get.temp.var.quantiles(filled.list$tmax, date_series, bs.date.series, temp.qtiles, bs.date.range, n, TRUE, min.base.data.fraction.present), assign.env = quantiles_env)
  }
  if ("tmin" %in% present_vars && length(intersect(tmin.dates, bs.date.series)) > days_threshold) {
    delayedAssign("tmin", climdex.pcic:::get.temp.var.quantiles(filled.list$tmin, date_series, bs.date.series, temp.qtiles, bs.date.range, n, TRUE, min.base.data.fraction.present), assign.env = quantiles_env)
  }
  if ("prec" %in% present_vars && length(intersect(prec.dates, bs.date.series)) > days_threshold) {
    delayedAssign("prec", climdex.pcic:::get.prec.var.quantiles(filled.list$prec, date_series, bs.date.range, prec.qtiles), assign.env = quantiles_env)
  }
  
  # --- 6. Construct the intermediate climdexInput object ---
  climdex_input_object <- new("climdexInput", data = filled.list, quantiles = quantiles_env, namasks = namasks,
                              dates = date_series, jdays = clim_data$yday, base.range = bs.date.range,
                              date.factors = date_factors, northern.hemisphere = northern.hemisphere,
                              max.missing.days = max.missing.days)
  
  #-------------------------------------------------------------------------------#
  # PART 2: CALCULATE INDICES AND RETURN FINAL RESULT
  #-------------------------------------------------------------------------------#
  
  # --- 7. Decide what to return ---
  if (!calculate_indices) {
    # If FALSE, return the intermediate object for other uses.
    return(climdex_input_object)
  }
  
  # --- 8. Calculate all ETCCDI indices ---
  ETCCDI <- list(
    FD      = function(x) climdex.fd(x),
    SU      = function(x) climdex.su(x),
    ID      = function(x) climdex.id(x),
    TR      = function(x) climdex.tr(x),
    GSL     = function(x) climdex.gsl(x, gsl.mode="GSL"),
    TXx     = function(x) climdex.txx(x, freq="annual"),
    TNx     = function(x) climdex.tnx(x, freq="annual"),
    TXn     = function(x) climdex.txn(x, freq="annual"),
    TNn     = function(x) climdex.tnn(x, freq="annual"),
    TN10p   = function(x) climdex.tn10p(x, freq="annual"),
    TX10p   = function(x) climdex.tx10p(x, freq="annual"),
    TN90p   = function(x) climdex.tn90p(x, freq="annual"),
    TX90p   = function(x) climdex.tx90p(x, freq="annual"),
    WSDI    = function(x) climdex.wsdi(x),
    CSDI    = function(x) climdex.csdi(x),
    DTR     = function(x) climdex.dtr(x, freq="annual"),
    Rx1day  = function(x) climdex.rx1day(x, freq="annual"),
    Rx5day  = function(x) climdex.rx5day(x, freq="annual"),
    SDII    = function(x) climdex.sdii(x),
    R5mm    = function(x) climdex.rnnmm(x, threshold=5),
    R10mm   = function(x) climdex.r10mm(x),
    R20mm   = function(x) climdex.r20mm(x),
    CDD     = function(x) climdex.cdd(x),
    CWD     = function(x) climdex.cwd(x),
    R95pTOT = function(x) climdex.r95ptot(x),
    R99pTOT = function(x) climdex.r99ptot(x),
    PRCPTOT = function(x) climdex.prcptot(x)
  )
  
  results_list <- lapply(ETCCDI, function(f) f(climdex_input_object))
  
  # Return the final list of calculated index results.
  return(results_list)
}

#=================================================================================#
# EXAMPLE USAGE
#=================================================================================#
# cat("--- Running Example ---\n")
#
# # 1. Generate sample data
# start_date <- as.Date("1981-01-01")
# end_date <- as.Date("2020-12-31")
# full_dates <- seq(start_date, end_date, by = "day")
# dates_no_leap <- full_dates[!(format(full_dates, "%m-%d") == "02-29")]
# pcict_dates_40yr <- PCICt::as.PCICt(as.character(dates_no_leap), cal = "365_day")
# n_days <- length(pcict_dates_40yr)
# tmx_40yr <- round(runif(n_days, min = 15, max = 35), 2)
# tmn_40yr <- tmx_40yr - round(runif(n_days, min = 5, max = 15), 2)
# pcp_40yr <- round(rexp(n_days, rate = 0.5), 2)
# base_period_40yr <- c(1981, 2010)
#
# # 2. Call the all-in-one function
# etccdi_results <- calculate_all_etccdi(
#   tmax       = tmx_40yr,
#   tmin       = tmn_40yr,
#   prec       = pcp_40yr,
#   tmax.dates = pcict_dates_40yr,
#   tmin.dates = pcict_dates_40yr,
#   prec.dates = pcict_dates_40yr,
#   base.range = base_period_40yr
# )
#
# # 3. View the results
# cat("Calculation complete. Viewing the structure of the results list:\n")
# str(etccdi_results, max.level = 1)