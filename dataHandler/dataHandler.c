#include "dataHandler.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

const double kb = 1.380649e-23;
const double pi = 3.14159265358979323846;
const double h_plank = 6.62607015e-34;
const double c = 299792458.0;
const double one_atm = 101325.0;

static unsigned int count_lines(char *content) {
    unsigned int count = 0;
    while (*content) {  // while the content pointer is on a valid char, continue, if we find the "\0" char, we stop
        char *eol = strchr(content, '\n');  // eol take the value of the pointer pointing on the first occurence of '\n' in the char[] content
        if (!eol) break;  // if eol is the null pointer, ther is no more '\n' in the file so we have parsed it all and stop the parsing
        if (eol - content > 10) count++;  // filters empty lines (the last line of the file for instance). (eol - content) returns the diff in position between the currnent reading head and the end of the line. If the line is big enough, we count it
        content = eol + 1;  // we set the reading head position just after the end of line to continue looping on lines
    }
    return count;
}


// converts a 36 base char into integer
// '0'-'9' → 0-9, 'A'-'Z' → 10-35, 'a'-'z' → 10-35
static int char_base36(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'Z') return c - 'A' + 10;
    if (c >= 'a' && c <= 'z') return c - 'a' + 10;
    return 0;
}

// converts a string of n char into 36 base
static int parse_base36(const char *s, int n) {
    int result = 0;
    for (int i = 0; i < n; i++)
        result = result * 36 + char_base36(s[i]);
    return result;
}

// reads n chars from src, converts into double. The inline keyword pasts the function code where it is called
static inline double parse_double(const char *src, int n) {
    char buf[32];
    memcpy(buf, src, n);  // copies the n first chars of src into buff
    buf[n] = '\0';  // terminates the line correctly
    return strtod(buf, NULL);  // converts char[] to double
}

static inline int parse_int(const char *src, int n) {
    char buf[16];
    memcpy(buf, src, n);  // copies the n first chars of src into buff
    buf[n] = '\0';  // terminates the line correctly
    return (int)strtol(buf, NULL, 10);  // converts char[] to int in base 10
}

static void readHITRANLine(char *line, HITRANStructure *data) {
    // manual parsing of the first 2 lines in base 36
    data->molecule_id    = parse_base36(line,     2);  // cols 1-2
    data->local_iso_id = parse_base36(line + 2, 1);  // col  3

    // continues manual parsing for speed
    data->nu = parse_double(line +  3, 12);
    data->sw = parse_double(line + 15, 10);
    data->a = parse_double(line + 25, 10);
    data->gamma_air = parse_double(line + 35,  5);
    data->gamma_self = parse_double(line + 40,  5);
    data->elower = parse_double(line + 45, 10);
    data->n_air = parse_double(line + 55,  4);
    data->delta_air = parse_double(line + 59,  8);
    memcpy(data->global_upper_quanta, line + 67, 15); data->global_upper_quanta[15] = '\0';
    memcpy(data->global_lower_quanta, line + 82, 15); data->global_lower_quanta[15] = '\0';
    memcpy(data->local_upper_quanta, line + 97, 15); data->local_upper_quanta[15]  = '\0';
    memcpy(data->local_lower_quanta, line + 112, 15); data->local_lower_quanta[15] = '\0';
    memcpy(data->ierr, line + 127,  6); data->ierr[6] = '\0';
    memcpy(data->iref, line + 133, 12); data->iref[12] = '\0';
    memcpy(data->line_mixing_flag, line + 145,  1); data->line_mixing_flag[1] = '\0';
    data->gp = parse_double(line + 146,  7);
    data->gpp = parse_double(line + 153,  7);
    // automatic parsing of the remaining (this function is kind of slow)
    //sscanf(line + 3, "%12lf%10lf%10lf%5lf%5lf%10lf%4lf%8lf%15c%15c%15c%15c%6c%12c%1c%7lf%7lf",
    //    &data->nu, &data->sw, &data->a,
    //    &data->gamma_air, &data->gamma_self,
    //    &data->elower, &data->n_air, &data->delta_air,
    //    data->global_upper_quanta, data->global_lower_quanta,
    //    data->local_upper_quanta,  data->local_lower_quanta,
    //    data->ierr, data->iref, data->line_mixing_flag,
    //    &data->gp, &data->gpp);  // parsing the line like it is prescribes by HITRAN
        // sscanf does not terminate the char variables properly (it only reads the text) so we must terminate the strings properly
    //    data->global_upper_quanta[15] = '\0';
    //    data->global_lower_quanta[15] = '\0';
    //    data->local_upper_quanta[15]  = '\0';
    //    data->local_lower_quanta[15]  = '\0';
    //    data->ierr[6]             = '\0';
    //    data->iref[12]            = '\0';
    //    data->line_mixing_flag[1] = '\0';

}

static void check_file_init(FILE *file) {
    // if the pointer *file == 0 (NULL pointer), it means that
    // the CURL object is not pointing to any valid location, so that it is not defined
    if (!file) {
        // printing in the standard error pipe (terminal by default) that the file init failed
        fprintf(stderr, "FILE init failed\n");
        // stopping the code with exit code 1
        exit(1);
    }
}

HITRANStructure *readHITRANData(char *path, unsigned int *array_size) {

    // opening the file to be parsed
    FILE *file = fopen(path, "r");
    check_file_init(file);

    // getting the number of char in the file
    fseek(file, 0, SEEK_END);  // puts the cursor at a 0 offset from the SEEK_END (so the end of tthe file). SO it puts the cursor at the end
    long file_size = ftell(file);

    // put the cursor back to the top of the file
    rewind(file);

    // reading the file in one function call instead of line by line to improve speed (but at the price of increased memory usage as the file must be stored entirely in memory)
    char *content = malloc(file_size + 1);  // allocating a char long enough to read the full file
    fread(content, 1, file_size, file);  // read the full file (fgets reads a line, fread reads the whole file)
    content[file_size] = '\0';  // terminating the char variable correctly
    fclose(file);  // closing the file

    // computing line number
    *array_size = count_lines(content);

    // allocating memory (allocating an array big enough to contain all the data structures)
    HITRANStructure *full_data = malloc(sizeof(HITRANStructure) * (*array_size));

    // we parse the file line by line
    char *p = content;
    unsigned int i = 0;
    p = content;
    while (*p) {  // read teh count_line function for details
        char *eol = strchr(p, '\n');
        if (!eol) break;
        if (eol - p > 10)
            readHITRANLine(p, &full_data[i++]);  // i is incremented AFTER the acces to full_data. If we write ++i, it is incremented before
        p = eol + 1;
    }


    // cleaning
    free(content);

    return(full_data);
}

HITRANStructure *getIsotopologByID(HITRANStructure *data, unsigned int *array_size, unsigned int iso_id) {
    //TODO: add a sorting algo juste to be sure that nu is ordered

    // counting the number of lines to keep
    unsigned int occ = 0;
    for (unsigned int i = 0; i < *array_size; i++) {
        if (data[i].local_iso_id == iso_id) {
            occ++;
        }
    }

    // Defining the result struct
    HITRANStructure *isotopologData = malloc(sizeof(HITRANStructure) * occ);  // cannot init with *isotopologData[occ] as it creates in the stack instead of the heap and the memory is freed at function end

    // copying the right lines into the new struct
    unsigned int j = 0;
    for (unsigned int i = 0; i < *array_size; i++) {
        if (data[i].local_iso_id == iso_id) {
            isotopologData[j] = data[i];  // isotopologData++ = data[i] allows only one loop instead of 2
            j++;
        }
    }

    *array_size = occ;

    return(isotopologData);
}

// the restrict keyword is here for automatic vectorization
double *getTau(double *nu, unsigned int len_nu, HITRANStructure *data, unsigned int len_data, double T, double P, double speConc) {

    // defining reference values
    double T0 = 296.0;
    double P_pa = P * one_atm;
    double n_mol = (P_pa / (kb * T)) * pow(1e-2, 3);

    // initianilze tau
    double *tau = calloc(len_nu, sizeof(double));  // calloc initialize to 0 while malloc does not initialize

    // computing gamma
    //double gamma[len_data]; // allocated like this, memory on the stack and can overflow if large data
    double *gamma = malloc(len_data * sizeof(double));  // allocated like this, on the heap, no stack overflow. It is longer to allocate but same access speed

    double T0_over_T = T0 / T;  // precomputing constants for slight speed improvement
    double P_times_1mc = P * (1.0 - speConc);  // precomputing constants for slight speed improvement
    double P_times_c   = P * speConc;  // precomputing constants for slight speed improvement

    for (unsigned int i = 0; i < len_data; i++) {
        gamma[i] = pow(T0_over_T, data[i].n_air) * (data[i].gamma_air * P_times_1mc + data[i].gamma_self * P_times_c);
    }

    // computing shifted line center
    //double nu_shifted[len_data];  // allocated like this, memory on the stack and can overflow if large data
    double *nu_shifted = malloc(len_data * sizeof(double));  // allocated like this, on the heap, no stack overflow. It is longer to allocate but same access speed


    for (unsigned int i = 0; i < len_data; i++) {
        nu_shifted[i] = data[i].nu + P * data[i].delta_air;
    }

    // computing Lorentzian profile (to improve computation speed, we can compute tau on a moving windowas contributions decrease in (delta nu) ^ 2
    for (unsigned int j = 0; j < len_data; j++) {  // looping on j first improves the cach efficiency by reading gamma value once for the whole i loop
        double nu_s = nu_shifted[j];
        double gam  = gamma[j];
        double sw_j = data[j].sw;
        double gam2 = gam * gam;
        double fL_cst = (1.0 / pi) * gam;

        for (unsigned int i = 0; i < len_nu; i++) {
            double delta_nu = nu[i] - nu_s;
            tau[i] += fL_cst / (gam2 + delta_nu * delta_nu) * sw_j;
        }
    }

    for (unsigned int i = 0; i < len_nu; i++) {
        tau[i] *= n_mol * speConc;
    }

    free(gamma);
    free(nu_shifted);
    return(tau);

}

double *plankProfile(double *nu, unsigned int len_nu, double T) {

    // initialize the profile
    double *profile = malloc(len_nu * sizeof(double));

    // computing exponential factor
    double *x = malloc(len_nu * sizeof(double));
    for (unsigned int i = 0; i < len_nu; i++) {
        x[i] = h_plank * c * 100.0 * nu[i] / (kb * T);
    }

    // computing the profile without overflow limiters
    for (unsigned int i = 0; i < len_nu; i++) {
        profile[i] = (2.0 * h_plank * pow(c, 2.0) * 1e8 * pow(nu[i], 3.0) / (exp(x[i]) - 1.0));
    }
    free(x);
    return(profile);

}

double *alpha(double *nu, unsigned int len_nu, HITRANStructure *data_CH4, unsigned int len_data_CH4,
              HITRANStructure *data_H2O, unsigned int len_data_H2O,
              HITRANStructure *data_CO2, unsigned int len_data_CO2,
              HITRANStructure *data_N2O, unsigned int len_data_N2O,
              AtmProp *layer_data) {

    // initialising alpha
    double *alph = calloc(len_nu, sizeof(double));

    // adding species by species contribution
    double *tau = getTau(nu, len_nu, data_CH4, len_data_CH4, layer_data->T, layer_data->P / one_atm, layer_data->C_CH4);
    for (unsigned int i = 0; i < len_nu; i++) {
        alph[i] += tau[i];
    }

    free(tau);

    tau = getTau(nu, len_nu, data_H2O, len_data_H2O, layer_data->T, layer_data->P / one_atm, layer_data->C_H2O);
    for (unsigned int i = 0; i < len_nu; i++) {
        alph[i] += tau[i];
    }
    free(tau);

    tau = getTau(nu, len_nu, data_N2O, len_data_N2O, layer_data->T, layer_data->P / one_atm, layer_data->C_N2O);
    for (unsigned int i = 0; i < len_nu; i++) {
        alph[i] += tau[i];
    }

    free(tau);

    tau = getTau(nu, len_nu, data_CO2, len_data_CO2, layer_data->T, layer_data->P / one_atm, layer_data->C_CO2);
    for (unsigned int i = 0; i < len_nu; i++) {
        alph[i] += tau[i];
    }

    free(tau);

    return(alph);

}




// Example of an approximate tau and alpha solver

double *getTau_approx(double *nu, unsigned int len_nu,
                      HITRANStructure *data, unsigned int len_data,
                      double T, double P, double speConc,
                      double tol) {

    const double T0          = 296.0;
    const double P_pa        = P * one_atm;
    const double n_mol       = (P_pa / (kb * T)) * pow(1e-2, 3);
    const double T0_over_T   = T0 / T;
    const double P_times_1mc = P * (1.0 - speConc);
    const double P_times_c   = P * speConc;

    double *tau        = calloc(len_nu, sizeof(double));
    double *gamma      = malloc(len_data * sizeof(double));
    double *nu_shifted = malloc(len_data * sizeof(double));

    for (unsigned int j = 0; j < len_data; j++) {
        gamma[j]      = pow(T0_over_T, data[j].n_air)
                      * (data[j].gamma_air * P_times_1mc + data[j].gamma_self * P_times_c);
        nu_shifted[j] = data[j].nu + P * data[j].delta_air;
    }

    for (unsigned int j = 0; j < len_data; j++) {

        const double gam    = gamma[j];
        const double gam2   = gam * gam;
        const double fL_cst = (1.0 / pi) * gam;
        const double sw_j   = data[j].sw;
        const double nu_s   = nu_shifted[j];
        const double ref    = fL_cst / gam2 * sw_j;

        // Binary search of the closest point to nu_s
        unsigned int lo = 0, hi = len_nu - 1;
        while (lo < hi) {
            unsigned int mid = (lo + hi) / 2;
            if (nu[mid] < nu_s) lo = mid + 1;
            else                hi = mid;
        }
        if (lo > 0 && fabs(nu[lo - 1] - nu_s) < fabs(nu[lo] - nu_s))
            lo--;
        unsigned int i0 = lo;

        // Center contribution
        double delta = nu[i0] - nu_s;
        tau[i0] += fL_cst / (gam2 + delta * delta) * sw_j;

        // Symetrical expansion
        int          left         = (int)i0 - 1;
        unsigned int right        = i0 + 1;
        int          left_active  = (left >= 0);
        int          right_active = (right < len_nu);

        while (left_active || right_active) {

            if (left_active) {
                delta = nu[left] - nu_s;
                double cc = fL_cst / (gam2 + delta * delta) * sw_j;
                tau[left] += cc;
                if (cc < tol * ref) left_active = 0;
                else { left--; if (left < 0) left_active = 0; }
            }

            if (right_active) {
                delta = nu[right] - nu_s;
                double cc = fL_cst / (gam2 + delta * delta) * sw_j;
                tau[right] += cc;
                if (cc < tol * ref) right_active = 0;
                else { right++; if (right >= len_nu) right_active = 0; }
            }
        }
    }

    for (unsigned int i = 0; i < len_nu; i++)
        tau[i] *= n_mol * speConc;

    free(gamma);
    free(nu_shifted);
    return tau;
}


double *alpha_approx(double *nu, unsigned int len_nu,
                     HITRANStructure *data_CH4, unsigned int len_data_CH4,
                     HITRANStructure *data_H2O, unsigned int len_data_H2O,
                     HITRANStructure *data_CO2, unsigned int len_data_CO2,
                     HITRANStructure *data_N2O, unsigned int len_data_N2O,
                     AtmProp *layer_data,
                     double tol) {

    double *alph = calloc(len_nu, sizeof(double));

    double *tau = getTau_approx(nu, len_nu, data_CH4, len_data_CH4, layer_data->T, layer_data->P / one_atm, layer_data->C_CH4, tol);
    for (unsigned int i = 0; i < len_nu; i++) alph[i] += tau[i];
    free(tau);

    tau = getTau_approx(nu, len_nu, data_H2O, len_data_H2O, layer_data->T, layer_data->P / one_atm, layer_data->C_H2O, tol);
    for (unsigned int i = 0; i < len_nu; i++) alph[i] += tau[i];
    free(tau);

    tau = getTau_approx(nu, len_nu, data_N2O, len_data_N2O, layer_data->T, layer_data->P / one_atm, layer_data->C_N2O, tol);
    for (unsigned int i = 0; i < len_nu; i++) alph[i] += tau[i];
    free(tau);

    tau = getTau_approx(nu, len_nu, data_CO2, len_data_CO2, layer_data->T, layer_data->P / one_atm, layer_data->C_CO2, tol);
    for (unsigned int i = 0; i < len_nu; i++) alph[i] += tau[i];
    free(tau);

    return alph;
}