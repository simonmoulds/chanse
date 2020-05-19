#include <RcppArmadillo.h>
#include <numeric>
#include "gwmodel.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// See here for how to document Rcpp function - http://r-pkgs.had.co.nz/src.html

// [[Rcpp::export]]
List model(NumericVector prec,	                 // precipitation, daily
	   double rain_rand,
           NumericVector reference_ET,           // reference ET (Hargreaves-Samani), daily
           IntegerVector day,	                 // Julian day, daily
           IntegerVector year,	                 // year, daily
	   IntegerVector month,	                 // month, daily
	   NumericVector Kc,                     // crop coefficient, daily
	   NumericVector et_depletion_factor,    // et depletion factor, daily
	   NumericVector rooting_depth,          // rooting depth, daily
	   NumericVector diesel_price,		 // diesel price, daily
	   NumericVector rice_price_per_t,       // price per tonne, daily
	   NumericVector wheat_price_per_t,      // price per tonne, daily
	   NumericVector rice_yield,             // tonne / hectare, annual
	   NumericVector wheat_yield,            // tonne / hectare, annual
           double rice_yield_coeff,
           double wheat_yield_coeff,
           IntegerMatrix wheat_irr_schedule,     // wheat irrigation schedule
           IntegerMatrix wheat_can_irr_schedule, // canal irrigation schedule for wheat
           IntegerMatrix rice_irr_schedule,	 // rice irrigation schedule
           IntegerMatrix rice_can_irr_schedule,  // canal irrigation schedule for rice
	   IntegerVector wheat_months,           // months of the year when wheat is grown
	   IntegerVector rice_months,            // months of the year when rice is grown
           double wheat_min,                     // minimum wheat irrigation depth
           double wheat_max,                     // maximum wheat irrigation depth
           double rice_min,                      // minimum rice irrigation depth
           double rice_max,                      // maximum rice irrigation depth
	   double Ky_wheat,			 // wheat yield response factor (FAO66)
	   double Ky_rice,			 // rice yield response factor (FAO66)
	   NumericMatrix fertiliser_price,       // price per kg, daily
	   NumericMatrix fertiliser_app_rate,    // fertiliser application rate kg / hectare, daily
	   LogicalVector farmer_canal_access,    // whether a farmer has access to a canal or not
	   IntegerVector init_farmer_category,   // initial farmer categories
           NumericVector well_depth,             // depth of each well
           double max_well_depth,                // maximum possible well depth
           NumericVector well_cost,              // cost of installing a new well - daily time series
           double initial_H,	                 // initial water table elevation TODO: multiply by Sy (unlike orig)
           double Sy,		                 // specific yield
           double initial_root_zone_depletion,   // initial root zone depletion
	   int harvest_day,                      
           int n_year,
           int n_day,
           int n_crop,
           int n_farmer,
           double region_area,                   // total area of all farms
           double farm_area,                     // area of each farm TODO: make this a vector, to allow different sizes
           double surface_elev,                  
	   double return_flow_coef,          
	   double irrigation_eff,               
	   double canal_volume,
	   double canal_leakage_coef,
	   double field_capacity,
	   double wilting_point,
	   double runoff_coef,
	   double Ks_max,
	   double Ks_min,
	   double evaporation_loss_coef,
           double fuel_eff,
           double saving_percentage,
           int saving_lower,
           int saving_upper) {

  double root_zone_depletion_init = initial_root_zone_depletion;
  double H_init = initial_H;
  NumericVector reference_crop_ET = reference_ET * Kc;
  
  // TODO: make crop_calendar a LogicalMatrix input argument, ncol=ncrop, nrow=n_day
  LogicalVector wheat_crop_cal(n_day);
  LogicalVector rice_crop_cal(n_day);
  for (int i = 0; i < n_day; i++) {
    wheat_crop_cal[i] = std::find(wheat_months.begin(), wheat_months.end(), month[i]) != wheat_months.end();
    rice_crop_cal[i] = std::find(rice_months.begin(), rice_months.end(), month[i]) != rice_months.end();
  }
  
  // get fertiliser price/app rate
  NumericVector fert_n_t_ha = fertiliser_app_rate(_,0);
  NumericVector fert_p_t_ha = fertiliser_app_rate(_,1);
  NumericVector fert_k_t_ha = fertiliser_app_rate(_,2);
  NumericVector fert_n_price_per_kg = fertiliser_price(_,0);
  NumericVector fert_p_price_per_kg = fertiliser_price(_,1);
  NumericVector fert_k_price_per_kg = fertiliser_price(_,2);

  // Preallocate output matrices
  IntegerMatrix farmer_category(n_year,n_farmer);
  // NumericMatrix farmer_category(n_year,n_farmer);
  NumericMatrix farmer_saving(n_year,n_farmer);
  NumericMatrix farmer_wheat_yield(n_year,n_farmer);
  NumericMatrix farmer_rice_yield(n_year,n_farmer);
  NumericMatrix farmer_wheat_irr_cost(n_year,n_farmer);
  NumericMatrix farmer_rice_irr_cost(n_year,n_farmer);
  NumericMatrix farmer_livelihood(n_year,n_farmer);
  NumericMatrix farmer_actual_ET(n_day,n_farmer);
  NumericMatrix farmer_gw_head(n_day,n_farmer);
  NumericMatrix farmer_root_zone_depletion(n_day,n_farmer);
  NumericMatrix farmer_abstraction(n_day,n_farmer);
  NumericMatrix farmer_outflow(n_day,n_farmer);
  NumericMatrix farmer_recharge(n_day,n_farmer);
  NumericMatrix farmer_irr_cost(n_day,n_farmer);

  // set initial values
  farmer_category(0,_) = init_farmer_category;
  // double H = H_init;
  // double root_zone_depletion = root_zone_depletion_init; 

  // start loop
  int count = 0;  
  for (int i = 0; i < n_day; i++) {
    for (int j = 0; j < n_farmer; j++) {
      
      // ================================
      // irrigation scheduling
      // ================================

      // Variables wd and rd represent the irrigation depth applied to wheat and
      // rice, respectively
      double wd = runif(1, wheat_min, wheat_max)[0];
      double rd = runif(1, rice_min, rice_max)[0];

      // the variable wheat|rice_irr_schedule is a matrix in which rows represent
      // the days on which irrigation takes place and columns represent individual
      // farmers. 
      IntegerVector wheat_irr_schedule_j = wheat_irr_schedule(_,j);
      IntegerVector rice_irr_schedule_j = rice_irr_schedule(_,j);

      // Variable wheat|rice_can_irr_schedule has the same structure as
      // wheat|rice_irr_schedule but indicates the days on which canal water is
      // available 
      IntegerVector wheat_can_irr_schedule_j = wheat_can_irr_schedule(_,j);
      IntegerVector rice_can_irr_schedule_j = rice_can_irr_schedule(_,j);

      // the purpose of this section is to find out whether irrigation takes
      // place on the current day. It does this by searching for the current day
      // in the irrigation schedule for the respective crops (i.e. wheat, rice)
      // and corresponding canal irrigation schedule.

      // wheat - total
      bool test = std::find(wheat_irr_schedule_j.begin(), wheat_irr_schedule_j.end(), day[i]) != wheat_irr_schedule_j.end();
      double irr_wh = 0;
      if (test) {
        irr_wh = wd;
      }

      // wheat - canal
      test = std::find(wheat_can_irr_schedule_j.begin(), wheat_can_irr_schedule_j.end(), day[i]) != wheat_can_irr_schedule_j.end();
      double irr_wh_canal = 0;
      if (test) {
        irr_wh_canal = wd;
      }

      // rice - total
      test = std::find(rice_irr_schedule_j.begin(), rice_irr_schedule_j.end(), day[i]) != rice_irr_schedule_j.end();
      double irr_ri = 0;
      if (test) {
        irr_ri = rd;
      }

      // rice - canal
      test = std::find(rice_can_irr_schedule_j.begin(), rice_can_irr_schedule_j.end(), day[i]) != rice_can_irr_schedule_j.end();
      double irr_ri_canal = 0;
      if (test) {
        irr_ri_canal = rd;
      }

      // TODO: rename these variables to something more informative
      double f1 = Rcpp::max(NumericVector::create(irr_wh,irr_ri)); // total irrigation
      double c1 = Rcpp::max(NumericVector::create(irr_wh_canal,irr_ri_canal)); // canal irrigation
      
      // ================================
      // groundwater model
      // ================================

      // the groundwater model runs every day and returns a list containing crop
      // ET, groundwater head and root zone depletion.
      List gw_sim = gwmodel(prec[i] * rain_rand,                  // precipitation on day i
			    reference_crop_ET[i],     // ET on day i
			    et_depletion_factor[i],   // et depletion factor on day i
			    rooting_depth[i],         // rooting depth on day i
			    H_init,                   // previous groundwater head
			    root_zone_depletion_init, // root zone depletion
			    f1,
			    c1,
			    well_depth[farmer_category(count,j) - 1],
			    max_well_depth,
			    farmer_canal_access[j],
			    farm_area,
			    surface_elev,
			    return_flow_coef,
			    irrigation_eff,
			    canal_volume,
			    canal_leakage_coef,
			    field_capacity,
			    wilting_point,
			    runoff_coef,
			    Ks_max,
			    Ks_min,
			    evaporation_loss_coef);

      farmer_actual_ET(i,j) = gw_sim["actual_ET"];
      farmer_gw_head(i,j) = gw_sim["H"];
      farmer_root_zone_depletion(i,j) = gw_sim["root_zone_depletion"];
      farmer_abstraction(i,j) = gw_sim["abstraction"];
      farmer_outflow(i,j) = gw_sim["outflow"];
      farmer_recharge(i,j) = gw_sim["recharge"];              
      // ================================
      // irrigation cost model
      // ================================

      // there needs to be an option for electricity
      
      // depth of groundwater from the surface (i.e. the depth from which water
      // needs to be pumped)
      double gw_depth = (surface_elev / region_area) - (farmer_gw_head(i,j) / region_area);

      double fuel = 0;
      if (c1 == 0) {
        fuel = (0.1133 * gw_depth) + 0.7949;
      }

      double fuel_m3 = fuel * fuel_eff / 102.87;              // fuel consumed to pump 1m3 water
      double fuel_ha = farmer_abstraction(i,j) * fuel_m3;     // volume of fuel used
      farmer_irr_cost(i,j) = fuel_ha * diesel_price[i] * 0.4; // irrigation cost per hectare

      // ================================
      // crop yield model
      // ================================

      if (day[i] == harvest_day) {

	// TODO: include labour factor
	
	LogicalVector index = (((year == (year[i] - 1)) & (day >= harvest_day)) | ((year == year[i]) & (day < harvest_day)));
	LogicalVector wheat_index = index & wheat_crop_cal;
	LogicalVector rice_index = index & rice_crop_cal;

	NumericVector farmer_actual_ET_j = farmer_actual_ET(_,j);
        NumericVector aet_wheat_daily = farmer_actual_ET_j[wheat_index];
        NumericVector ref_et_wheat_daily = reference_crop_ET[wheat_index];
	double aet_wheat = sum(aet_wheat_daily);
	double ref_et_wheat = sum(ref_et_wheat_daily);

        NumericVector aet_rice_daily = farmer_actual_ET_j[rice_index];
        NumericVector ref_et_rice_daily = reference_crop_ET[rice_index];
	double aet_rice = sum(aet_rice_daily);
	double ref_et_rice = sum(ref_et_rice_daily);

	// QUERY: what does (rice|wheat)_yield represent - maximum yield?
	if (aet_wheat > 0) {
	  farmer_wheat_yield(count,j) = (1 - (Ky_wheat * (1 - (aet_wheat / ref_et_wheat)))) * wheat_yield[count] * wheat_yield_coeff;
	} else {
	  farmer_wheat_yield(count,j) = 0;
	}

	if (aet_rice > 0) {
	  farmer_rice_yield(count,j) = (1 - (Ky_rice * (1 - (aet_rice / ref_et_rice)))) * rice_yield[count] * rice_yield_coeff;
	} else {
	  farmer_rice_yield(count,j) = 0;
	}

        NumericVector farmer_irr_cost_j = farmer_irr_cost(_,j);
        NumericVector wheat_irr_cost_daily = farmer_irr_cost_j[wheat_index];
        NumericVector rice_irr_cost_daily = farmer_irr_cost_j[rice_index];
        farmer_wheat_irr_cost(count,j) = sum(wheat_irr_cost_daily);
        farmer_rice_irr_cost(count,j) = sum(rice_irr_cost_daily);

        // get the price of rice and wheat at harvest time, assuming that
        // harvest is the last day of the crop season

        // NB JoK was taking the mean price over the season - check
        // whether this is right
        
        IntegerVector v = seq(0, wheat_index.size() - 1);
        IntegerVector vv = v[wheat_index];
        int wheat_index_last = tail(vv, 1)[0];

        v = seq(0, rice_index.size() - 1);
        vv = v[rice_index];
        int rice_index_last = tail(vv, 1)[0];
          
        double wheat_harvest_price_per_t = wheat_price_per_t[wheat_index_last];
        double rice_harvest_price_per_t = rice_price_per_t[rice_index_last];

        // get mean fertiliser application rate over the growing season
        NumericVector n_app_rate = fert_n_t_ha[index];
        NumericVector p_app_rate = fert_p_t_ha[index];
        NumericVector k_app_rate = fert_k_t_ha[index];
        
        double mean_fert_n_t_ha = mean(n_app_rate);
        double mean_fert_p_t_ha = mean(p_app_rate);
        double mean_fert_k_t_ha = mean(k_app_rate);

        NumericVector n_price = fert_n_price_per_kg[index];
        NumericVector p_price = fert_p_price_per_kg[index];
        NumericVector k_price = fert_k_price_per_kg[index];
        
        double mean_fert_n_price_per_kg = mean(n_price);
        double mean_fert_p_price_per_kg = mean(p_price);
        double mean_fert_k_price_per_kg = mean(k_price);

        // assume we double because there are two growing seasons?
        double total_fert_cost = (mean_fert_n_t_ha * mean_fert_n_price_per_kg * farm_area / 10000 +
                                  mean_fert_p_t_ha * mean_fert_p_price_per_kg * farm_area / 10000 +
                                  mean_fert_k_t_ha * mean_fert_k_price_per_kg * farm_area / 10000) * 2;

        double total_irrigation_cost = farmer_wheat_irr_cost(count,j) + farmer_rice_irr_cost(count,j);
	
        double total_income = (farmer_wheat_yield(count,j) * wheat_harvest_price_per_t + farmer_rice_yield(count,j) * rice_harvest_price_per_t);
	
        farmer_livelihood(count,j) = total_income - total_irrigation_cost - total_fert_cost;

	// savings for 1st run calculated here (as a % of income) ~~~ MAKE SURE SAVINGS MATCH THAT IN TOP SECTION
	// if (count == 0) {
        //   double rand = round(runif(1, saving_lower, saving_upper)[0]);
	//   farmer_saving(count,j) = farmer_livelihood(count,j) * saving_percentage + rand;
	// } else {
	//   farmer_saving(count,j) = farmer_livelihood(count,j) * saving_percentage + farmer_saving(count-1,j);

	//   // check (i) farmer can afford new well, (ii) farmer can change category
	//   bool cond = ((farmer_saving(count,j) > well_cost[i]) & (farmer_category(count-1,j) < 3));

	//   if (cond) {
	//     farmer_category(count,j) = farmer_category(count-1,j) + 1;
	//     farmer_saving(count,j) = farmer_saving(count,j) - well_cost[i];
	//   } else {
	//     farmer_category(count,j) = farmer_category(count-1,j);
	//   }
	// }
	if (count == 0) {
          double rand = round(runif(1, saving_lower, saving_upper)[0]);
	  farmer_saving(count,j) = farmer_livelihood(count,j) * saving_percentage + rand;
	} else {
	  farmer_saving(count,j) = farmer_livelihood(count,j) * saving_percentage + farmer_saving(count-1,j);
        }

        // check (i) farmer can afford new well, (ii) farmer can change category
        bool cond = ((farmer_saving(count,j) > well_cost[i]) && (farmer_category(count,j) < 3));
        // bool cond = ((farmer_saving(count,j) > well_cost[i]) && (farmer_category(count-1,j) < 3));

        if (cond) {
          farmer_saving(count,j) = farmer_saving(count,j) - well_cost[i];
	}

	if ((count+1) < n_year) {
 	  if (cond) {
	    farmer_category(count+1,j) = farmer_category(count,j) + 1;
	  } else {
            farmer_category(count+1,j) = farmer_category(count,j); // ERROR - something to do with 'count+1'
	  }
        }
      }

      // update groundwater head & rzd for next farmer
      // NB this approach implies that the farmers access the water in order, so
      // farmer 1 will access shallower water than farmer 2, and so on.
      H_init = farmer_gw_head(i,j); // TODO: change this - all farmers should have the same initial value (?)
      root_zone_depletion_init = farmer_root_zone_depletion(i,j); // TODO: change this - rzd should be particular to individual farmer
      
    }

    // update count at harvest day
    if ((day[i] == harvest_day) & ((count+1) < n_year)) {
      count++;
    }
  }

  return List::create(Named("farmer_category") = farmer_category,
                      _["farmer_livelihood"] = farmer_livelihood,
                      _["farmer_saving"] = farmer_saving,
                      _["farmer_wheat_yield"] = farmer_wheat_yield,
                      _["farmer_rice_yield"] = farmer_rice_yield,
                      _["farmer_actual_ET"] = farmer_actual_ET,
                      _["farmer_gw_head"] = farmer_gw_head,
                      _["farmer_root_zone_depletion"] = farmer_root_zone_depletion,
                      _["farmer_abstraction"] = farmer_abstraction,
                      _["farmer_outflow"] = farmer_outflow,
                      _["farmer_recharge"] = farmer_recharge);
  
}
