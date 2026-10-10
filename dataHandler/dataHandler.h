#pragma once
#include "../athmosphericLayers/US_standard_atm/standardAtm.h"

typedef struct {
    /*
     * Data structure containing the HITRAN data
     */
    unsigned int molecule_id;
    unsigned int local_iso_id;
    double nu;
    double sw;
    double a;
    double gamma_air;
    double gamma_self;
    double elower;
    double n_air;
    double delta_air;
    double gp;
    double gpp;
    char   global_upper_quanta[16];
    char   global_lower_quanta[16];
    char   local_upper_quanta[16];
    char   local_lower_quanta[16];
    char   ierr[7];
    char   iref[13];
    char   line_mixing_flag[2];
} HITRANStructure;


HITRANStructure *readHITRANData(char *path, unsigned int *array_size);  // parsing the whole file and returning the appropriate structure
HITRANStructure *getIsotopologByID(HITRANStructure *data, unsigned int *array_size, unsigned int iso_id);

double *getTau(double *nu, unsigned int len_nu, HITRANStructure *data, unsigned int len_data, double T, double P, double speConc);

double *plankProfile(double *nu, unsigned int len_nu, double T);

double *alpha(double *nu, unsigned int len_nu, HITRANStructure *data_CH4, unsigned int len_data_CH4,
              HITRANStructure *data_H2O, unsigned int len_data_H2O,
              HITRANStructure *data_CO2, unsigned int len_data_CO2,
              HITRANStructure *data_N2O, unsigned int len_data_N2O,
              AtmProp *layer_data);

// Example of an approximate tau and alpha solver

double *getTau_approx(double *nu, unsigned int len_nu, HITRANStructure *data, unsigned int len_data,
                      double T, double P, double speConc, double tol);

double *alpha_approx(double *nu, unsigned int len_nu,
                     HITRANStructure *data_CH4, unsigned int len_data_CH4,
                     HITRANStructure *data_H2O, unsigned int len_data_H2O,
                     HITRANStructure *data_CO2, unsigned int len_data_CO2,
                     HITRANStructure *data_N2O, unsigned int len_data_N2O,
                     AtmProp *layer_data, double tol);