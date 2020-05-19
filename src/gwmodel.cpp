#include <RcppArmadillo.h>
#include <numeric>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// See here for how to document Rcpp function - http://r-pkgs.had.co.nz/src.html

// [[Rcpp::export]]
List gwmodel(double prec,                     // precipitation (mm)
	     double reference_crop_ET,        // reference crop evapotranspiration (mm)
	     double et_depletion_factor,      // ET depletion factor (FAO56)
	     double rooting_depth,	      // rooting depth (mm)
	     double H_init,                   // groundwater head from previous time point
	     double root_zone_depletion_init, // root zone depletion
	     double f1,		              // total irrigation [demand?]
	     double c1,		              // canal irrigation [demand?]
	     double well_depth,               // current well depth (m)
	     double max_well_depth,           // maximum well depth (m)
	     bool farmer_canal_access,        // access to canal water?
	     double farm_area,                // farm area (m2)
	     double surface_elev,             // surface elevation
	     double return_flow_coef,         // return flow coefficient (-)
	     double irrigation_eff,           // irrigation efficiency (-)
	     double canal_volume,             // canal volume, depth * width * length (m3)
	     double canal_leakage_coef,       // canal leakage (-) 
	     double field_capacity,           // field capacity, (mm)
	     double wilting_point,            // wilting point (mm)
	     double runoff_coef,              // runoff coefficient (-)
	     double Ks_max,                   // max water stress coefficient (-)
	     double Ks_min,                   // min water stress coefficient (-)
	     double evaporation_loss_coef) {

  // infiltration to the root zone
  double infiltration_volume = (prec / 1000) * runoff_coef * evaporation_loss_coef * farm_area;

  double irrigation_volume;
  if ((farmer_canal_access) || (H_init > well_depth)) {
    irrigation_volume = (farm_area * f1) * (1 + irrigation_eff); // * evaporation_loss_coef;
  } else {
    irrigation_volume = 0;
  }

  // QUERY: is return flow a separate mechanism from recharge?
  // double return_flow_volume = irrigation_volume * return_flow_coef;
  // double return_flow_volume = 0;

  // total and readily available water content
  double taw = 1000 * ((field_capacity - wilting_point) * rooting_depth); // m -> mm
  double raw = et_depletion_factor * taw;

  // crop stress, actual ET, root zone depletion
  double Ks;
  if (root_zone_depletion_init < raw) {
    Ks = 1.0;
  } else {
    Ks = (taw - root_zone_depletion_init) / ((1 - et_depletion_factor) * taw);
  }

  double actual_ET = Ks * reference_crop_ET;
  
  // NB abstraction converted to depth (by dividing by farm_area) and
  // units changed to mm
  double root_zone_depletion_new = root_zone_depletion_init + actual_ET - (((irrigation_volume * evaporation_loss_coef + infiltration_volume) / farm_area) * 1000);
  
  double recharge_depth = 0.0;
  if (root_zone_depletion_new < 0) {
    recharge_depth = (root_zone_depletion_new / 1000) * -1; // * baseflow_index ???
    root_zone_depletion_new = 0;
  }
  double recharge_volume = recharge_depth * farm_area;
    
  // hydraulic head after removing drainage due to lateral flow                                                                     
   double k = 30; // 30; // 0.086;
   double z = 79 - 10;
   double L = 11500;
   double w = 100;
   double b = (79 - (H_init/10000)) - (79 - 30);  
   double cond = ((k * w * b) /L) /10000;
   double outflow =((((cond * (H_init/10000 - z)))*-1));
   double H_new_lat = (H_init - outflow);
  
 // double H_new_lat = cond*(H_init - z);                                                                                        //  double H_new_lat = H_init;                                                                                                 // Rprintf("H_init: %f\n", H_init);  
  // Rprintf("H_new_lat: %f\n", H_new_lat);                                                                                     //  Rprintf("outflow  : %f\n", outflow);                                                                                      //   Rprintf("b  : %f\n", b);
  // Rprintf("cond     : %f\n", cond);
  // Rprintf("H_init   : %f/n", H_init);    
  // Rprintf("surface_elev   : %f/n", surface_elev);
  // Rprintf("k   : %f/n", k);
  // Rprintf("b   : %f/n",b);

  // hydraulic head
  double H_new;
  double canal_leakage = 0;
  if (farmer_canal_access) {

    // NB concerned about canal leakage - why divide by farm_area??? old line below: H_new = H_init + recharge_volume + (canal_volume * canal_leakage_coef / farm_area);    
    canal_leakage = (canal_volume * canal_leakage_coef);
    // canal_leakage = (canal_volume * canal_leakage_coef / farm_area);
    // Rprintf("canal leakage : %f\n", canal_leakage);
    // Rprintf("recharge vol  : %f\n", recharge_volume);      
    H_new = H_new_lat + recharge_volume + canal_leakage;

  } else {

    if (H_new_lat <= well_depth) {
      
      // in this case the farmer's well is not deep enough to exploit the resource,
      // so the new groundwater head is simply that for the previous day plus
     // recharge
      H_new = H_new_lat + recharge_volume;

    } else {

      // otherwise the farmer's well is deep enough to extract groundwater
      H_new = H_new_lat + recharge_volume - irrigation_volume;
    }
  }

  // update recharge volume to include canal leakage
  recharge_volume += canal_leakage;
  
  // groundwater level cannot be above the ground surface
  if (H_new > surface_elev) {
    H_new = surface_elev; 
  }

  // nor, for the purposes of the model, can it be below the maximum well depth
  if (H_new <= max_well_depth) {
    H_new = max_well_depth;
  }
  
  return List::create(Named("actual_ET") = actual_ET, _["H"] = H_new, _["root_zone_depletion"] = root_zone_depletion_new, _["abstraction"] = irrigation_volume, _["outflow"] = outflow, _["recharge"] = recharge_volume);

}
