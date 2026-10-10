#include "../dataHandler/dataHandler.h"
#include "../athmosphericLayers/US_standard_atm/standardAtm.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>  // creating folder and having system infos

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
    // variable for selecting the specie we want to see
    const char *spe = "CO2";  // sould not be changed to prevent issues in alpha computation
    const unsigned int iso_id = 1;
    const unsigned int len_nu = 1000;
    const float nu_min = 0.0;
    const float nu_max = 10000.0;
    const float z = 0.0;

    // making the path in which to look for the data
    char path[2048];
    snprintf(path, sizeof(path), "data/bandes_%s.par", spe);

    // reading the data file and creating the HITRANStructure objcet
    unsigned int bandes_nb_CO2;
    HITRANStructure *bandes_CO2 = readHITRANData(path, &bandes_nb_CO2);

    // extracting isotopolog
    unsigned int istop_bandes_nb_CO2 = bandes_nb_CO2;
    HITRANStructure *solo_CO2 = getIsotopologByID(bandes_CO2, &istop_bandes_nb_CO2, iso_id);

    // printing the number of lines to be sure the inplace modif works right
    printf("bandes_nb = %d\n", bandes_nb_CO2);
    printf("istop_bandes_nb = %d\n", istop_bandes_nb_CO2);

    // getting tau
    double nu[len_nu];
    for (int i = 0; i < len_nu; i++) {
        nu[i] = nu_min + (nu_max - nu_min) * i / (len_nu - 1);
    }
    double *tau = getTau(nu,len_nu, solo_CO2, istop_bandes_nb_CO2, 300.0, 1.0, 420e-6);

    // getting plank profile
    double *plank = plankProfile(nu, len_nu, 300.0);

    // getting alpha
    AtmProp atmProp;
    get_atm_prop(z, &atmProp);
    unsigned int bandes_nb_H2O;
    HITRANStructure *bandes_H2O = readHITRANData("data/bandes_H2O.par", &bandes_nb_H2O);
    unsigned int istop_bandes_nb_H2O = bandes_nb_H2O;
    HITRANStructure *solo_H2O = getIsotopologByID(bandes_H2O, &istop_bandes_nb_H2O, iso_id);

    unsigned int bandes_nb_CH4;
    HITRANStructure *bandes_CH4 = readHITRANData("data/bandes_CH4.par", &bandes_nb_CH4);
    unsigned int istop_bandes_nb_CH4 = bandes_nb_CH4;
    HITRANStructure *solo_CH4 = getIsotopologByID(bandes_CH4, &istop_bandes_nb_CH4, iso_id);

    unsigned int bandes_nb_N2O;
    HITRANStructure *bandes_N2O = readHITRANData("data/bandes_N2O.par", &bandes_nb_N2O);
    unsigned int istop_bandes_nb_N2O = bandes_nb_N2O;
    HITRANStructure *solo_N2O = getIsotopologByID(bandes_N2O, &istop_bandes_nb_N2O, iso_id);

    double *alph = alpha(nu,len_nu, solo_CH4, istop_bandes_nb_CH4, solo_H2O, istop_bandes_nb_H2O, solo_CO2, istop_bandes_nb_CO2, solo_N2O, istop_bandes_nb_N2O, &atmProp);

    // definition of the output file for data visualization
    // creating saving folder if it does not exist
    #ifdef _WIN32  // on a windows computer
        mkdir("test_data_handler");
    #else  // on a linux computer
        mkdir("test_data_handler", 0777);  // giving the right access
    #endif

    FILE *file = fopen("test_data_handler/tau.txt", "w");
    check_file_init(file);
    for (int i = 0; i < len_nu; i++) {
        fprintf(file, "%e,%e,%e,%e\n", nu[i], tau[i], plank[i], alph[i]);
    }
    fclose(file);

    // launching the python plotting program
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "python \"%s/tests/test_data_handler.py\"", PROJECT_SOURCE_DIR);
    system(cmd);

    return 0;
}