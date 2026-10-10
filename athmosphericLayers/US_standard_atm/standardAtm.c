#include <math.h>  // using math functions
#include <stdlib.h>  // using exit()
#include <stdio.h>  // for prints
#include "standardAtm.h"


////////////////////
// atm definition //
////////////////////

const double T_STD = 288.15;  // mean see level temperature (K)
const double P_SDT = 101325.0;  // mean see level pressure (Pa)
const double RHO_STD = 1.125;  // mean see level density (kg/m3)
const double g = 9.80665;  // gravitational constant (m/s2)
const double r = 287.05;  // specific gaz constant for air (J/kg.K)
const double Re = 6356766.0;  // Mean Earth radius (m)

const AtmosphericLayer USstandardAtm[] = {
    /*
     * Definition of US standard atmosphere
     */
    {0.0, 11000.0, -0.0065, 288.15, 101325.0},  // Troposphere
    {11000.0, 20000.0, 0.0000, 216.65, 22632.0},  // Stratosphere 1
    {20000.0, 32000.0, 0.0010, 216.65, 5474.9},  // Stratosphere 2
    {32000.0, 47000.0, 0.0028, 228.65, 868.02},  // Stratosphere 3
    {47000.0, 51000.0, 0.0000, 270.65, 110.91},  // Stratosphere 4
    {51000.0, 71000.0, -0.0028, 270.65, 66.94},  // Mésosphere 1
    {71000.0, 86000.0, -0.0020, 214.65, 3.96}  // Mésosphere 2
};
const unsigned int nStandardAtmLayers = sizeof(USstandardAtm) / sizeof(AtmosphericLayer);

///////////////////////
// Utility functions //
///////////////////////

double h2z(double h) {
    /*
     * Converts geopotential height to geometric altitude
     * :param h: Geopotential height in meter
     * :return: Geometric altitude in meter
     */
    return Re * h / (Re - h);
}

double z2h(double z) {
    /*
     * Converts geometric altitude to geopotential height
     * :param z: Geometric altitude in meter
     * :return: Geopotential altitude in meter
     */
    return Re * z / (Re + z);
}

/////////////////////
// species profile //
/////////////////////

static double CO2(double h) {
    /*
     * Model for CO2 concentration. Based on constant concentration profile
     * :param h: Geopotential altitude in meter
     * :return: CO2 mol fraction
     */
    return 420e-6;
}

static double CH4(double h) {
    /*
     * Model for CH4 concentration. Based on constant concentration profile
     * :param h: Geopotential altitude in meter
     * :return: CH4 mol fraction
     */
    return 1.9e-6;
}

static double N2O(double h) {
    /*
     * Model for N2O concentration. Based on constant concentration profile
     * :param h: Geopotential altitude in meter
     * :return: N2O mol fraction
     */
    return 0.33e-6;
}

static double H2O(double h) {
    /*
     * Model for H2O concentration. Based on constant concentration profile
     * :param h: Geopotential altitude in meter
     * :return: H2O mol fraction
     */
    double C0 = 0.004;  // ground concentration of H2O (mol)
    double H = 2.0e3;  // characteristic length (m)
    return C0 * exp(- h / H);
}

///////////////////
// main function //
///////////////////

void get_atm_prop(double z, AtmProp *out) {
    /*
     * Computes atmosphere properties given geometric altitude
     * :param z: Geometric altitude in meter
     * :return out: in place modification of a struct containing Temperature, pressure and density of the atmosphere at geometric altitude computed with US Standard Atmosphere model
     */
    double h = z2h(z);  // converting geometric to geopotential

    // checking that the resulting geopotetial height is valid for the model
    if (h < 0.0){fprintf(stderr, "Impossible to compute standard atmosphere values for negative geopotential altitudes (h = %f)", h); exit(1);}
    if (h > 86000.0){fprintf(stderr, "Impossible to compute standard atmosphere values above a geopotential altitude of 86000 m (h = %f) with the US standard atmosphere model", h); exit(1);}

    for (unsigned int i = 0; i < nStandardAtmLayers; i++) {  // looping through each layer
        if (USstandardAtm[i].h_min <= h & h < USstandardAtm[i].h_max) {  // when we find the layer in which the computation altitude belongs, we compute atm props and store it in out

            double dh = h - USstandardAtm[i].h_min;  // we compute the altitude in the specific layer

            out->T = USstandardAtm[i].T_base + USstandardAtm[i].T_gradient * dh;  // we compute and store T

            if (USstandardAtm[i].T_gradient == 0) {  // if the layer is isothermal, the pressure has a specific formula
                out->P = USstandardAtm[i].P_base * exp(- g * dh / (r * USstandardAtm[i].T_base));  // we compute the pressure
            }
            else {  // if the layer is not isothermal, the pressure has a specific formula
                out->P = USstandardAtm[i].P_base * pow((out->T / USstandardAtm[i].T_base), - g / (USstandardAtm[i].T_gradient * r));
            }

            out->rho = out->P / (r * out->T);

            out->C_CH4 = CH4(h);    // computing species concentration
            out->C_H2O = H2O(h);
            out->C_CO2 = CO2(h);
            out->C_N2O = N2O(h);

            return; // stop the function as there is no mean in continuing the search

        }
    }

}