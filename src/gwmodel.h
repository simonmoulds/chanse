#ifndef GWMODEL_H
#define GWMODEL_H

#include <Rcpp.h>
Rcpp::List gwmodel(double prec,
		   double reference_crop_ET,
		   double et_depletion_factor,
		   double rooting_depth,
		   double H_init,
		   double root_zone_depletion_init,
		   double f1,
		   double c1,
		   double well_depth,
		   double max_well_depth,
		   bool farmer_canal_access,
		   double farm_area,
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
		   double evaporation_loss_coef);
#endif
