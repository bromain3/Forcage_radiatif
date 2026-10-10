#pragma once
#include "../athmosphericLayers/US_standard_atm/standardAtm.h"
#include "../dataHandler/dataHandler.h"
#include <math.h>

void propagate_through_layer_inplace(double *I, unsigned int len_I,
                                     double *nu, unsigned int len_nu,
                                     HITRANStructure *data_CH4, unsigned int len_data_CH4,
                                     HITRANStructure *data_H2O, unsigned int len_data_H2O,
                                     HITRANStructure *data_CO2, unsigned int len_data_CO2,
                                     HITRANStructure *data_N2O, unsigned int len_data_N2O,
                                     double z_base, double dz);

double *propagate_inplace(double *nu, unsigned int len_nu,double z_min, double z_max, unsigned int n_z,
                          HITRANStructure *data_CH4, unsigned int len_data_CH4,
                          HITRANStructure *data_H2O, unsigned int len_data_H2O,
                          HITRANStructure *data_CO2, unsigned int len_data_CO2,
                          HITRANStructure *data_N2O, unsigned int len_data_N2O);

void propagate_through_layer_inplace_approx(double *I, unsigned int len_I,
                                     double *nu, unsigned int len_nu,
                                     HITRANStructure *data_CH4, unsigned int len_data_CH4,
                                     HITRANStructure *data_H2O, unsigned int len_data_H2O,
                                     HITRANStructure *data_CO2, unsigned int len_data_CO2,
                                     HITRANStructure *data_N2O, unsigned int len_data_N2O,
                                     double z_base, double dz, double tol);

double *propagate_inplace_approx(double *nu, unsigned int len_nu,double z_min, double z_max, unsigned int n_z,
                          HITRANStructure *data_CH4, unsigned int len_data_CH4,
                          HITRANStructure *data_H2O, unsigned int len_data_H2O,
                          HITRANStructure *data_CO2, unsigned int len_data_CO2,
                          HITRANStructure *data_N2O, unsigned int len_data_N2O, double tol);