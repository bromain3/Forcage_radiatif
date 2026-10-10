#pragma once

extern const double T_STD;  // mean see level temperature (K)
extern const double P_SDT;  // mean see level pressure (Pa)
extern const double RHO_STD;  // mean see level density (kg/m3)
extern const double g;  // gravitational constant (m/s2)
extern const double r;  // specific gaz constant for air (J/kg.K)
extern const double Re;  // Mean Earth radius (m)

typedef struct {  //TODO: refactor this file in a "general" folder ? so that each atmospheric model can rely on this struct
    /*
     * Data structure defining all variable describing an atmospheric layer.
     */
    const double h_min;  // Minimal geopotential height of the layer in meter
    const double h_max;  // Maximal geopotential height of the layer in meter
    const double T_gradient;  // Temperature gradient of the layer in Kelvin/kilometer
    const double T_base;  // Temperature of the layer base in Kelvin
    const double P_base;  // Pressure of the layer base in Pa
} AtmosphericLayer;

typedef struct {
    /*
     * Data structure for storing all atmosphere properties
     */
    double T;
    double P;
    double rho;
    double C_CO2;
    double C_CH4;
    double C_N2O;
    double C_H2O;
} AtmProp;

void get_atm_prop(double z, AtmProp *out);

double h2z(double h);
double z2h(double h);