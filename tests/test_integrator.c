#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>  // creating folder and having system infos
#include "../integrator/integrator.h"

void check_file_init(FILE *file) {
    // if the pointer *file == 0 (NULL pointer), it means that
    // the CURL object is not pointing to any valid location, so that it is not defined
    if (!file) {
        // printing in the standard error pipe (terminal by default) that the file init failed
        fprintf(stderr, "FILE init failed\n");
        // stopping the code with exit code 1
        exit(1);
    }
}

int main() {

    // definition of parameters
    const unsigned int len_nu = 1000;
    const float nu_min = 200.0;
    const float nu_max = 1600.0;
    const float z = 0.0;
    const float dz = 15000.0;

    const float z_min = 0.0;
    const float z_max = 15000.0;
    const unsigned int n_z = 100;

    const double tol = 1.0e-6;


    // initialization of atm properties at z
    AtmProp atm_prop;
    get_atm_prop(z, &atm_prop);

    // initialization of nu list
    double nu[len_nu];
    for (int i = 0; i < len_nu; i++) {
        nu[i] = nu_min + (nu_max - nu_min) * i / (len_nu - 1);
    }

    // getting all optical data
    unsigned int bandes_nb_H2O;
    HITRANStructure *bandes_H2O = readHITRANData("data/bandes_H2O.par", &bandes_nb_H2O);
    unsigned int istop_bandes_nb_H2O = bandes_nb_H2O;
    HITRANStructure *solo_H2O = getIsotopologByID(bandes_H2O, &istop_bandes_nb_H2O, 1);

    unsigned int bandes_nb_CH4;
    HITRANStructure *bandes_CH4 = readHITRANData("data/bandes_CH4.par", &bandes_nb_CH4);
    unsigned int istop_bandes_nb_CH4 = bandes_nb_CH4;
    HITRANStructure *solo_CH4 = getIsotopologByID(bandes_CH4, &istop_bandes_nb_CH4, 1);

    unsigned int bandes_nb_N2O;
    HITRANStructure *bandes_N2O = readHITRANData("data/bandes_N2O.par", &bandes_nb_N2O);
    unsigned int istop_bandes_nb_N2O = bandes_nb_N2O;
    HITRANStructure *solo_N2O = getIsotopologByID(bandes_N2O, &istop_bandes_nb_N2O, 1);

    unsigned int bandes_nb_CO2;
    HITRANStructure *bandes_CO2 = readHITRANData("data/bandes_CO2.par", &bandes_nb_CO2);
    unsigned int istop_bandes_nb_CO2 = bandes_nb_CO2;
    HITRANStructure *solo_CO2 = getIsotopologByID(bandes_CO2, &istop_bandes_nb_CO2, 1);

    // computing blackbody radiation at alt z (used as initial conditions)
    double *I = plankProfile(nu, len_nu, atm_prop.T);

    // storing initial profile in a file
    // creating saving folder if it does not exist
    #ifdef _WIN32  // on a windows computer
        mkdir("test_integrator");
    #else  // on a linux computer
        mkdir("test_integrator", 0777);  // giving the right access
    #endif
    FILE *file = fopen("test_integrator/data.txt", "w");
    check_file_init(file);
    for (int i = 0; i < len_nu-1; i++) {  // writing nu
        fprintf(file, "%e,", nu[i]);
    }fprintf(file, "%e\n", nu[len_nu-1]);
    for (int i = 0; i < len_nu-1; i++) {  //writing I
        fprintf(file, "%e,", I[i]);
    }fprintf(file, "%e\n", I[len_nu-1]);

    // propagating radiation through the layer
    propagate_through_layer_inplace(I, len_nu,
                                     nu, len_nu,
                                     solo_CH4, istop_bandes_nb_CH4,
                                     solo_H2O, istop_bandes_nb_H2O,
                                     solo_CO2, istop_bandes_nb_CO2,
                                     solo_N2O, istop_bandes_nb_N2O,
                                     z, dz);

    for (int i = 0; i < len_nu-1; i++) {
        fprintf(file, "%e,", I[i]);
    }fprintf(file, "%e\n", I[len_nu-1]);

    // memory free of I for preventing memory leaks
    free(I);

    // testing the full script
    I = propagate_inplace(nu, len_nu, z_min, z_max, n_z,
                         solo_CH4, istop_bandes_nb_CH4,
                         solo_H2O, istop_bandes_nb_H2O,
                         solo_CO2, istop_bandes_nb_CO2,
                         solo_N2O, istop_bandes_nb_N2O);

    for (int i = 0; i < len_nu-1; i++) {
        fprintf(file, "%e,", I[i]);
    }fprintf(file, "%e", I[len_nu-1]);
    fclose(file);

    // launching the python plotting program
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "python \"%s/tests/test_integrator.py\"", PROJECT_SOURCE_DIR);
    system(cmd);

}