#include "integrator.h"

#include <stdio.h>
#include <stdlib.h>

void propagate_through_layer_inplace(double *I, unsigned int len_I,
                                     double *nu, unsigned int len_nu,
                                     HITRANStructure *data_CH4, unsigned int len_data_CH4,
                                     HITRANStructure *data_H2O, unsigned int len_data_H2O,
                                     HITRANStructure *data_CO2, unsigned int len_data_CO2,
                                     HITRANStructure *data_N2O, unsigned int len_data_N2O,
                                     double z_base, double dz) {

    // computing the point at which to evaluate optical properties (we use a midpoint scheme)
    double z_mid = z_base + dz / 2.0;
    double dz_cm = dz * 100.0;

    // computing atm properties at this layer
    AtmProp atm_prop;
    get_atm_prop(z_mid, &atm_prop);

    // computing alpha
    double *alph = alpha(nu,len_nu, data_CH4, len_data_CH4, data_H2O, len_data_H2O, data_CO2, len_data_CO2, data_N2O, len_data_N2O, &atm_prop);

    // computing black body radiation
    double *B = plankProfile(nu, len_nu, atm_prop.T);

    // in place modification of I
    for (unsigned int i = 0; i < len_I; i++) {
        double e_tau = exp(-alph[i] * dz_cm);
        I[i] = I[i] * e_tau + B[i] * (1.0 - e_tau);
    }

    free(alph);
    free(B);

}


double *propagate_inplace(double *nu, unsigned int len_nu,double z_min, double z_max, unsigned int n_z,
                          HITRANStructure *data_CH4, unsigned int len_data_CH4,
                          HITRANStructure *data_H2O, unsigned int len_data_H2O,
                          HITRANStructure *data_CO2, unsigned int len_data_CO2,
                          HITRANStructure *data_N2O, unsigned int len_data_N2O) {

    // getting atm properties at z_min (boundary condition)
    AtmProp atm_prop;
    get_atm_prop(z_min, &atm_prop);

    // initializing irradiance at z_min
    double *I = plankProfile(nu, len_nu, atm_prop.T);

    // computing dz (static step)
    double dz = (z_max - z_min) / (double)n_z;
    double z = z_min;
    printf("Atmosphere divided in %d steps (%.2f m per step)\n", n_z, dz);

    for (unsigned int i = 0; i < n_z; i++) {

        if (i % (n_z / 20) == 0) {
            printf("\tz = %.2f m  (%.2f %%)\n", z, 100.0 * (double)i / (double)n_z);
        }

        propagate_through_layer_inplace(I, len_nu,
                                 nu, len_nu,
                                 data_CH4, len_data_CH4,
                                 data_H2O, len_data_H2O,
                                 data_CO2, len_data_CO2,
                                 data_N2O, len_data_N2O,
                                 z, dz);
        z += dz;

    }

    return(I);

}




// Example of an approximate tau and alpha solver

void propagate_through_layer_inplace_approx(double *I, unsigned int len_I,
                                     double *nu, unsigned int len_nu,
                                     HITRANStructure *data_CH4, unsigned int len_data_CH4,
                                     HITRANStructure *data_H2O, unsigned int len_data_H2O,
                                     HITRANStructure *data_CO2, unsigned int len_data_CO2,
                                     HITRANStructure *data_N2O, unsigned int len_data_N2O,
                                     double z_base, double dz, double tol) {

    // computing the point at which to evaluate optical properties (we use a midpoint scheme)
    double z_mid = z_base + dz / 2.0;
    double dz_cm = dz * 100.0;

    // computing atm properties at this layer
    AtmProp atm_prop;
    get_atm_prop(z_mid, &atm_prop);

    // computing alpha
    double *alph = alpha_approx(nu,len_nu, data_CH4, len_data_CH4, data_H2O, len_data_H2O, data_CO2, len_data_CO2, data_N2O, len_data_N2O, &atm_prop, tol);

    // computing black body radiation
    double *B = plankProfile(nu, len_nu, atm_prop.T);

    // in place modification of I
    for (unsigned int i = 0; i < len_I; i++) {
        double e_tau = exp(-alph[i] * dz_cm);
        I[i] = I[i] * e_tau + B[i] * (1.0 - e_tau);
    }

    free(alph);
    free(B);

}


double *propagate_inplace_approx(double *nu, unsigned int len_nu,double z_min, double z_max, unsigned int n_z,
                          HITRANStructure *data_CH4, unsigned int len_data_CH4,
                          HITRANStructure *data_H2O, unsigned int len_data_H2O,
                          HITRANStructure *data_CO2, unsigned int len_data_CO2,
                          HITRANStructure *data_N2O, unsigned int len_data_N2O, double tol) {

    // getting atm properties at z_min (boundary condition)
    AtmProp atm_prop;
    get_atm_prop(z_min, &atm_prop);

    // initializing irradiance at z_min
    double *I = plankProfile(nu, len_nu, atm_prop.T);

    // computing dz (static step)
    double dz = (z_max - z_min) / (double)n_z;
    double z = z_min;
    printf("Atmosphere divided in %d steps (%.2f m per step)\n", n_z, dz);

    for (unsigned int i = 0; i < n_z; i++) {

        if (i % (n_z / 20) == 0) {
            printf("\tz = %.2f m  (%.2f %%)\n", z, 100.0 * (double)i / (double)n_z);
        }

        propagate_through_layer_inplace_approx(I, len_nu,
                                 nu, len_nu,
                                 data_CH4, len_data_CH4,
                                 data_H2O, len_data_H2O,
                                 data_CO2, len_data_CO2,
                                 data_N2O, len_data_N2O,
                                 z, dz, tol);
        z += dz;

    }

    return(I);

}