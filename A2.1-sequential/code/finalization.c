#include <stdio.h>
#include "util_write_files.h"

void finalization(char* file_in, int nprocs, int myrank, int total_iters, double residual_ratio,
                  int nintci, int nintcf, double* var) {

    char file_out[100];
    sprintf(file_out, "%s_summary.r%d.out", file_in, myrank);

    int status = store_simulation_stats(file_in, file_out, nintci, nintcf, var, total_iters,
                                        residual_ratio);

    if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
}

