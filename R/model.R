#' Create a model object.
#'
#' @param weather_data weather_data. See \code{\link{weather_data}}.
#' @param price_data price_data. See \code{\link{price_data}}.
#' @param yield_data yield_data. See \code{\link{yield_data}}.
#' @param farm_parameters farm_parameters. See \code{\link{farm_parameters}}
#' @param crop_parameters crop_parameters. See \code{\link{crop_parameters}}
#' @param irrigation_schedule irrigation_schedule. See \code{\link{irrigation_schedule}}
#' @param initial_condition initial_condition. See \code{\link{initial_condition}}
#'
#' @return model
#' 
model = function(weather_data,
                 price_data,
                 yield_data,
                 farm_parameters,
                 crop_parameters,
                 irrigation_schedule,
                 initial_condition) {

    ## TODO: checks
    
    x = list(
        weather_data = weather_data,
        price_data = price_data,
        yield_data = yield_data,
        farm_parameters = farm_parameters,
        crop_parameters = crop_parameters,
        irrigation_schedule = irrigation_schedule,
        initial_condition = initial_condition
    )
    structure(x, "model")
}

## ## Author : Simon Moulds & Jimmy O'Keeffe
## ## Date   : January 2019

## library(Rcpp)
## library(magrittr)
## library(xts)
## sourceCpp("code/model.cpp")

## ## parameters
## ## ##########
## source(file.path("code", "config", params_fn))

## ## irrigation schedule
## ## ###################
## source(file.path("code", "config", irrigation_schedule_fn))

## ## region characteristics
## ## ######################

## ## total area of model region
## region_area = farm_area * n_farmer

## ## surface elevation
## surface_elev = surface_elev_masl * region_area * Sy

## ## initial head
## H_init = H_init_masl * region_area * Sy

## ## well depths of three categories
## well_depth1 = well_depth1_masl * region_area * Sy
## well_depth2 = well_depth2_masl * region_area * Sy
## well_depth3 = well_depth3_masl * region_area * Sy
## well_depth = c(well_depth1, well_depth2, well_depth3)    
## max_well_depth = min(well_depth)                         

## ## climate
## ## #######

## climate_data = read.table(file.path("data","input",climate_data_fn), header=TRUE)
## rainfall_data = read.table(file.path("data","input",rainfall_data_fn), header=TRUE)

## prec = rainfall_data[["precip"]]
## reference_ET = climate_data[["ETo_mm"]]

## ## time data - taken from time attributes of climate input data
## tm = as.POSIXlt(rainfall_data[["date"]]) # safe to use this in future runs as well?
## ## day = tm[["yday"]] + 1
## ## month = tm[["mon"]] + 1
## ## year = tm[["year"]] + 1900
## day = tm$yday + 1
## month = tm$mon + 1
## year = tm$year + 1900

## ## Socioeconomic data
## ## ##################

## price_data = read.table(file.path("data","input",price_data_fn), header = TRUE)
## price_data_year = as.POSIXlt(price_data[["year"]])$year + 1900
## ## price_data_year = price_data[["year"]] %>% as.POSIXlt %>% `[[`("year") %>% `+`(1900)
## index1 = match(sort(unique(year)), price_data_year) ## this index removes years outside the model period
## index2 = match(year, price_data_year[index1])       ## this index maps annual values to daily

## rice_price_per_t = price_data[["rice_p_per_t"]][index1][index2]
## wheat_price_per_t = price_data[["wheat_p_per_t"]][index1][index2]

## fert_n_t_ha = price_data[["fert_n_t_ha"]][index1][index2]
## fert_p_t_ha = price_data[["fert_p_t_ha"]][index1][index2]
## fert_k_t_ha = price_data[["fert_k_t_ha"]][index1][index2]

## fert_n_price_per_kg = price_data[["fert_n_p_per_kg"]][index1][index2]
## fert_p_price_per_kg = price_data[["fert_p_p_per_kg"]][index1][index2]
## fert_k_price_per_kg = price_data[["fert_k_p_per_kg"]][index1][index2]

## ## make fertiliser price/application rate a daily time series
## fertiliser_price = matrix(data=c(fert_n_price_per_kg, fert_p_price_per_kg, fert_k_price_per_kg), ncol=3)
## fertiliser_app_rate = matrix(data=c(fert_n_t_ha, fert_p_t_ha, fert_k_t_ha), ncol=3)

## diesel_price_data = read.table(file.path("data","input",diesel_price_data_fn), header = TRUE)
## diesel_price_data_year = diesel_price_data[["year"]] %>% as.POSIXlt
## index = diesel_price_data_year %in% tm
## diesel_price = diesel_price_data[["dp"]][index]

## well_cost_data = read.table(file.path("data","input",well_cost_data_fn), header = TRUE)
## well_cost_data_year = well_cost_data[["year"]] %>% as.POSIXlt
## index = well_cost_data_year %in% tm
## well_cost = well_cost_data[["well_cost"]][index]

## ## yield data
## ## ##########

## yield_data = read.table(file.path("data","input",yield_data_fn), header=TRUE)
## yield_data_year = yield_data[["Year"]]
## index = yield_data_year %in% year
## wheat_yield = yield_data[["wheat_th"]][index]
## rice_yield = yield_data[["rice_th"]][index]

## ## crop data 
## ## #########

## assign_daily_vals = function(..., day, default_val) {
##     ## Function to assign values depending on the Julian day of the year
##     ## (e.g. Kc, ET depletion factor, rooting depth)
##     dots=list(...)
##     x = rep(default_val, length(day))
##     for (i in 1:length(dots)) {
##         start = dots[[i]][1]
##         end = dots[[i]][2]
##         val = dots[[i]][3]
##         if (end < start) {
##             x[day >= start] = val
##             x[day <= end] = val
##         } else {
##             x[day >= start & day <= end] = val
##         }
##     }
##     x
## }

## ## crop coefficient
## Kc = assign_daily_vals(Kc1w, Kc2w, Kc3w, Kc4w, Kc1r, Kc2r, Kc3r, Kc4r, day=day, default_val=0)

## ## ET depletion factor
## et_depletion_factor = assign_daily_vals(et_depletion_factor_wheat, et_depletion_factor_rice, day=day, default_val=0)

## ## rooting depth
## rooting_depth = assign_daily_vals(rooting_depth_wheat, rooting_depth_rice, day=day, default_val=0.01)

## ## TODO: is this the best way of doing this?
## et_depletion_factor[et_depletion_factor < et_depletion_factor_fallow] = et_depletion_factor_fallow;
## rooting_depth[rooting_depth < rooting_depth_fallow] = rooting_depth_fallow;

## ## run model
## ## #########

## n_year = length(unique(year))
## n_day = length(day)
## n_crop = 2
## n_wheat_irr = length(wheat_irr_schedule_mu)
## n_rice_irr = length(rice_irr_schedule_sd)

## for (k in 1:n_runs) {

##     wheat_irr_schedule =
##         rnorm(n = n_wheat_irr * n_farmer, mean=wheat_irr_schedule_mu, sd=wheat_irr_schedule_sd) %>%
##         round %>%
##         `[<-`(. < 1, 1) %>%
##         matrix(nrow=n_wheat_irr, ncol=n_farmer)

##     rice_irr_schedule =
##         rnorm(n = n_rice_irr * n_farmer, mean=rice_irr_schedule_mu, sd=rice_irr_schedule_sd) %>%
##         round %>%
##         `[<-`(. < 1, 1) %>%
##         matrix(nrow=n_rice_irr, ncol=n_farmer)

##     n_wheat_can_irr = length(wheat_can_irr_schedule_mu)
##     n_rice_can_irr = length(rice_can_irr_schedule_mu)

##     wheat_can_irr_schedule =
##         rnorm(n = n_wheat_can_irr * n_farmer, mean=wheat_can_irr_schedule_mu, sd=wheat_can_irr_schedule_sd) %>%
##         round %>%
##         `[<-`(. < 1, 1) %>%
##         matrix(nrow=n_wheat_can_irr, ncol=n_farmer)

##     rice_can_irr_schedule =
##         rnorm(n = n_rice_can_irr * n_farmer, mean=rice_can_irr_schedule_mu, sd=rice_can_irr_schedule_sd) %>%
##         round %>%
##         `[<-`(. < 1, 1) %>%
##         matrix(nrow=n_rice_can_irr, ncol=n_farmer)

##     ## this function changes the irrigation schedule if an irrigation is
##     ## scheduled a short period (say 1 or 2 days) after canal water is
##     ## available. This is based on the assumption that if a farmer is due
##     ## to irrigation on Tuesday (for example), but in fact canal water is
##     ## available on the Monday, he will irrigated a day early in order to
##     ## make use of canal water (need to check what JoK thinks of this
##     ## approach)
##     for (i in 1:nrow(wheat_irr_schedule)) {
##         for (j in 1:ncol(wheat_irr_schedule)) {
##             v = wheat_irr_schedule[i,j]
##             diff = v - wheat_can_irr_schedule[,j]
##             ix = diff > 0 & diff < 1
##             if (any(ix)) {
##                 v = wheat_can_irr_schedule[ix,j]
##             }
##             wheat_irr_schedule[i,j] = v
##         }
##     }

##     init_farmer_category = sample(c(1,2,3), size=n_farmer, prob=cat_prob, replace=TRUE)

##     ## farmer_canal_access = sample(x=0:1, size=n_farmer, replace=TRUE, prob=c(1-canal_prob, canal_prob)) %>% as.logical
##     ## canal_prob = rnorm(1, canal_prob_mu, canal_prob_sd)

##     canal_prob = runif(1, canal_prob_lower, canal_prob_upper)
##     farmer_canal_access_ix = sample(seq_len(n_farmer), size=floor(n_farmer * canal_prob), replace=FALSE)
##     farmer_canal_access = rep(0, length.out=n_farmer)
##     farmer_canal_access[farmer_canal_access_ix] = 1

##     rain_rand = runif(1, 0.95, 1.05);

##     x = model(prec, rain_rand, reference_ET, day, year, month, Kc, et_depletion_factor, rooting_depth, diesel_price, rice_price_per_t, wheat_price_per_t, rice_yield, wheat_yield, rice_yield_coeff, wheat_yield_coeff, wheat_irr_schedule, wheat_can_irr_schedule, rice_irr_schedule, rice_can_irr_schedule, wheat_months, rice_months, wheat_min, wheat_max, rice_min, rice_max, Ky_wheat, Ky_rice, fertiliser_price, fertiliser_app_rate, farmer_canal_access, init_farmer_category, well_depth, max_well_depth, well_cost, H_init, Sy, root_zone_depletion_init, harvest_day, n_year, n_day, n_crop, n_farmer, region_area, farm_area, surface_elev, return_flow_coef, irrigation_eff, canal_volume, canal_leakage_coef, field_capacity, wilting_point, runoff_coef, Ks_max, Ks_min, ## k, z, L, w, b, cond, 
##     evaporation_loss_coef, fuel_eff, saving_percentage, saving_lower, saving_upper)

##     ## change the format of the output to comply with JoK original
##     ## data analysis scripts
##     farmer_output = array(data=NA, dim=c(n_year, n_farmer, 6)) 
##     farmer_output[,,1] = x[["farmer_category"]]
##     farmer_output[,,2] = x[["farmer_livelihood"]]
##     farmer_output[,,3] = x[["farmer_saving"]]
##     farmer_output[,,4] = NA
##     farmer_output[,,5] = x[["farmer_wheat_yield"]]
##     farmer_output[,,6] = x[["farmer_rice_yield"]]
##     ## farmer_output[,,7] = x[["farmer_outflow"]]
##     ## farmer_output[,,8] = x[["farmer_recharge"]]
##     save(farmer_output, file=file.path("data","output","RData",paste0(run_id, "_farmer_run_", k, ".RData")))

##     gw_output = array(data=NA, dim=c(n_day, n_farmer, 8))
##     gw_output[,,1] = x[["farmer_actual_ET"]]
##     gw_output[,,2] = x[["farmer_gw_head"]]
##     gw_output[,,3] = x[["farmer_root_zone_depletion"]]
##     gw_output[,,4] = NA
##     gw_output[,,5] = NA
##     gw_output[,,6] = x[["farmer_abstraction"]]
##     gw_output[,,7] = x[["farmer_outflow"]]
##     gw_output[,,8] = x[["farmer_recharge"]]
##     save(gw_output, file=file.path("data","output","RData",paste0(run_id, "_gw_run_", k, ".RData")))
## }

## ## ======================================
## ## write txt file output
## ## ======================================

## for (i in 1:n_runs){
    
##     ## farmer
##     filein <- file.path("data","output","RData",paste0(run_id, "_farmer_run_",i,".RData"))
##     farmer.in <- get(load(filein))
    
##     if(i == 1){
##         cat <- array(NA,dim=c(n_year,n_farmer,n_runs))
##         income <- array(NA,dim=c(n_year,n_farmer,n_runs))
##         savings <- array(NA,dim=c(n_year,n_farmer,n_runs))
##         yieldW <- array(NA,dim=c(n_year,n_farmer,n_runs))
##         yieldR <- array(NA,dim=c(n_year,n_farmer,n_runs))
##     }
    
##     cat[,,i] <- farmer.in[,,1]
##     income[,,i] <- farmer.in[,,2]
##     savings[,,i] <- farmer.in[,,3]
##     yieldW[,,i] <- farmer.in[,,5]
##     yieldR[,,i] <- farmer.in[,,6]
    
##     ## groundwater
##     filein <- file.path("data","output","RData",paste0(run_id, "_gw_run_", i, ".RData"))
##     gw.in <- get(load(filein))
##     if (i == 1){
##         gw <- array(NA,dim=c(n_day, n_farmer, n_runs))
##         abstraction <- array(NA,dim=c(n_day, n_farmer, n_runs))
##         outflow <- array(NA,dim=c(n_day, n_farmer, n_runs))
##         recharge <- array(NA,dim=c(n_day, n_farmer, n_runs))
##     }
    
##     gw[,,i] <- (gw.in[,,2] / region_area) / Sy
##     abstraction[,,i] <- gw.in[,,6]
##     outflow[,,i] <- gw.in[,,7]
##     recharge[,,i] <- gw.in[,,8]    
## }

## ## average over farmers
## gw = apply(gw, c(1,3), mean)

## abstraction = apply(abstraction, c(1,3), mean)
## outflow = apply(outflow, c(1,3), mean)
## recharge = apply(recharge, c(1,3), mean)

## ## do something here

## ## idea:
## ## 1. make objects into zoo|xts
## ## 2. compute sum (monthly)
## ## 3. compute sum (annual)
## abstraction_xts = xts(as.data.frame(abstraction), order.by=tm)
## abstraction_monthly = apply.monthly(abstraction_xts, FUN=colSums)
## abstraction_annual = apply.yearly(abstraction_xts, FUN=colSums)

## recharge_xts = xts(as.data.frame(recharge), order.by=tm)
## recharge_monthly = apply.monthly(recharge_xts, FUN=colSums)
## recharge_annual = apply.yearly(recharge_xts, FUN=colSums)

## outflow_xts = xts(as.data.frame(outflow), order.by=tm)
## outflow_monthly = apply.monthly(outflow_xts, FUN=colSums)
## outflow_annual = apply.yearly(outflow_xts, FUN=colSums)

## total_in_monthly = recharge_monthly
## total_out_monthly = outflow_monthly + abstraction_monthly
## deltas_monthly = total_in_monthly - total_out_monthly


## total_in = recharge_annual
## total_out = outflow_annual + abstraction_annual
## deltas = total_in - total_out

## rdr = (recharge_monthly) - (abstraction_monthly + outflow_monthly)   ## was / now minus -
## ardr = apply.yearly(rdr, FUN=colSums)

## income = apply(income, c(1,3), mean)
## yieldW = apply(yieldW, c(1,3), mean)
## yieldR = apply(yieldR, c(1,3), mean)

## gw_mean <- apply(gw, 1, mean)  ## mean of Gw for plotting
## gw_max  <- apply(gw, 1, FUN=function(x) quantile(x, 0.99))
## gw_min  <- apply(gw, 1, FUN=function(x) quantile(x, 0.01))
## gw_med  <- apply(gw, 1, median)

## deltas_mean <- apply(deltas, 1, mean)  ## mean of Deltas for plotting
## deltas_max  <- apply(deltas, 1, FUN=function(x) quantile(x, 0.99))
## deltas_min  <- apply(deltas, 1, FUN=function(x) quantile(x, 0.01))
## deltas_med  <- apply(deltas, 1, median)


## deltas_monthly_mean <- apply(deltas_monthly, 1, mean)  ## mean of monthly Deltas for plotting


## rdr_mean <- apply(rdr, 1, mean)  ## mean of Rdr for plotting
## rdr_max  <- apply(rdr, 1, FUN=function(x) quantile(x, 0.99))
## rdr_min  <- apply(rdr, 1, FUN=function(x) quantile(x, 0.01))
## rdr_med  <- apply(rdr, 1, median)

## ardr_mean <- apply(ardr, 1, mean)  ## mean of Ardr for plotting
## ardr_max  <- apply(ardr, 1, FUN=function(x) quantile(x, 0.99))
## ardr_min  <- apply(ardr, 1, FUN=function(x) quantile(x, 0.01))
## ardr_med  <- apply(ardr, 1, median)

## income_mean <- apply(income, 1, mean) 
## income_max  <- apply(income, 1, FUN=function(x) quantile(x, 0.99)) 
## income_min  <- apply(income, 1, FUN=function(x) quantile(x, 0.01))
## income_med  <- apply(income, 1, median) 

## yieldW_mean <- apply(yieldW, 1, mean) 
## yieldW_max <- apply(yieldW, 1, FUN=function(x) quantile(x, 0.99)) 
## yieldW_min <- apply(yieldW, 1, FUN=function(x) quantile(x, 0.01))
## yieldW_med <- apply(yieldW, 1, median)

## yieldR_mean <- apply(yieldR, 1, mean) 
## yieldR_max <- apply(yieldR, 1, FUN=function(x) quantile(x, 0.99)) 
## yieldR_min <- apply(yieldR, 1, FUN=function(x) quantile(x, 0.01))
## yieldR_med <- apply(yieldR, 1, median)

## ## gw_mean <- apply(gw[,1,], 1, mean)  ## mean of Gw for plotting
## ## gw_min  <- apply(gw[,1,], 1, min)   ## min of Gw for plotting
## ## gw_max  <- apply(gw[,1,], 1, max)   ## max of Gw for plotting
## ## gw_med  <- apply(gw[,1,], 1, median)   ## max of Gw for plotting

## ## income_mean <- apply(income[,n_farmer,], 1, mean) 
## ## income_max  <- apply(income[,n_farmer,], 1, max) 
## ## income_min  <- apply(income[,n_farmer,], 1, min) 
## ## income_med  <- apply(income[,n_farmer,], 1, median) 

## ## yieldW_mean <- apply(yieldW[,n_farmer,], 1, mean) 
## ## yieldW_max <- apply(yieldW[,n_farmer,], 1, max) 
## ## yieldW_min <- apply(yieldW[,n_farmer,], 1, min)
## ## yieldW_med <- apply(yieldW[,n_farmer,], 1, median)

## ## yieldR_mean <- apply(yieldR[,n_farmer,], 1, mean) 
## ## yieldR_max <- apply(yieldR[,n_farmer,], 1, max) 
## ## yieldR_min <- apply(yieldR[,n_farmer,], 1, min)
## ## yieldR_med <- apply(yieldR[,n_farmer,], 1, median)

## write.table(gw_min, file.path("data","output","txt",paste0("gw_min_", run_id, ".txt")))
## write.table(gw_max, file.path("data","output","txt",paste0("gw_max_", run_id, ".txt")))
## write.table(gw_med, file.path("data","output","txt",paste0("gw_med_", run_id, ".txt")))

## write.table(rdr_min, file.path("data","output","txt",paste0("rdr_min_", run_id, ".txt")))
## write.table(rdr_max, file.path("data","output","txt",paste0("rdr_max_", run_id, ".txt")))
## write.table(rdr_med, file.path("data","output","txt",paste0("rdr_med_", run_id, ".txt")))
## write.table(rdr_mean, file.path("data","output","txt",paste0("rdr_mean_", run_id, ".txt")))

## write.table(ardr_min, file.path("data","output","txt",paste0("ardr_min_", run_id, ".txt")))
## write.table(ardr_max, file.path("data","output","txt",paste0("ardr_max_", run_id, ".txt")))
## write.table(ardr_med, file.path("data","output","txt",paste0("ardr_med_", run_id, ".txt")))
## write.table(ardr_mean, file.path("data","output","txt",paste0("ardr_mean_", run_id, ".txt")))

## write.table(deltas_mean, file.path("data","output","txt",paste0("deltas_mean_", run_id, ".txt")))
## write.table(deltas_monthly_mean, file.path("data","output","txt",paste0("deltas_monthly_mean_", run_id, ".txt")))


## write.table(income_min, file.path("data","output","txt",paste0("income_min_", run_id, ".txt")))
## write.table(income_max, file.path("data","output","txt",paste0("income_max_", run_id, ".txt")))
## write.table(income_med, file.path("data","output","txt",paste0("income_med_", run_id, ".txt")))

## write.table(yieldW_min, file.path("data","output","txt",paste0("yieldW_min_", run_id, ".txt")))
## write.table(yieldW_max, file.path("data","output","txt",paste0("yieldW_max_", run_id, ".txt")))
## write.table(yieldW_med, file.path("data","output","txt",paste0("yieldW_med_", run_id, ".txt")))

## write.table(yieldR_min, file.path("data","output","txt",paste0("yieldR_min_", run_id, ".txt")))
## write.table(yieldR_max, file.path("data","output","txt",paste0("yieldR_max_", run_id, ".txt")))
## write.table(yieldR_med, file.path("data","output","txt",paste0("yieldR_med_", run_id, ".txt")))
