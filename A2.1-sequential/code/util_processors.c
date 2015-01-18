#include "util_processors.h"
#include "posl_definitions.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


//TODO:use this function instead of any other command line arguments processor (like in gccg.c)
void process_cl(char* file_in, char* part_type, char* read_type, int *input_key, int *part_key, int *read_key) {
    if (strstr(file_in, "tjunc") != NULL) {
        *input_key = POSL_DATA_TJUNC;
    } else if (strstr(file_in, "drall") != NULL) {
        *input_key = POSL_DATA_DRALL;
    } else if (strstr(file_in, "pent") != NULL) {
        *input_key = POSL_DATA_PENT;
    } else if (strstr(file_in, "cojack") != NULL) {
        *input_key = POSL_DATA_COJACK;
    }

    if (!strcmp(part_type, "classic")) {
        *part_key = POSL_PARTITIONING_CLASSIC;
    } else if (!strcmp(part_type, "dual")) {
        *part_key = POSL_PARTITIONING_DUAL;
    } else {
        *part_key = POSL_PARTITIONING_NODAL;
    }

    if (!strcmp(read_type, "oneread")) {
        *read_key = POSL_INIT_ONE_READ;
    } else {
        *read_key = POSL_INIT_ALL_READ;
    }
}