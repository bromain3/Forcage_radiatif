#include <stdio.h>  // prints
#include <stdlib.h>  // using exit()
#include <sys/stat.h>  // creating folder and having system infos
#include "../athmosphericLayers/US_standard_atm/standardAtm.h"

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

void write_headers(FILE *file) {
    fputs("z,T,P,rho,CH4,CO2,N2O,H2O\n",file);
}

void write_line(FILE *file, char write_buffer[]) {
    fputs(write_buffer, file);
}

int main() {
    // definition of a structure to store atmosphere values
    AtmProp atmProp;

    // definition of number and range of atmosphere layers
    double z_min = 0.0;
    double z_max = 86000.0;
    double z;
    double steps = 1000.0;

    // define var for printing / saving
    char write_buffer[512];  // defines max length of data

    // definition of the output file for data visualization
    // creating saving folder if it does not exist
    #ifdef _WIN32  // on a windows computer
        mkdir("test_US_atm");
    #else  // on a linux computer
        mkdir("test_US_atm", 0777);  // giving the right access
    #endif
    FILE *file = fopen("test_US_atm/data.txt", "w");
    check_file_init(file);
    write_headers(file);

    for (unsigned int i = 0; i <= steps; i++) {

        z = z_min + (z_max - z_min) * i / steps;

        get_atm_prop(z, &atmProp);

        snprintf(write_buffer, sizeof(write_buffer), "%e,%e,%e,%e,%e,%e,%e,%e\n", z, atmProp.T, atmProp.P, atmProp.rho, atmProp.C_CH4, atmProp.C_CO2, atmProp.C_N2O, atmProp.C_H2O);

        write_line(file, write_buffer);

        //printf("Altitude: %6.0f m | Temperature: %6.2f K | Pression: %10.2f Pa | Densite: %.5f kg/m3\n", z, atmProp.T, atmProp.P, atmProp.rho);

    }
    // closing the file so that other programms can read it
    fclose(file);

    // launching the python plotting program
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "python \"%s/tests/test_US_std_atm.py\"", PROJECT_SOURCE_DIR);
    system(cmd);

    return 0;
}