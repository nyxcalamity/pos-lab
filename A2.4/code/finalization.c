/**
 * Finalization step - write results and other computational vectors to files
 *
 * @author V. Petkov, Denys Sobchyshak, Denys Korzh
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "util_write_files.h"
#include "posl_definitions.h"

void finalization(char* file_in, int nprocs, int myrank, int total_iters, double residual_ratio,
                  int nintci, int nintcf, double* var, int* local_global_index, int* global_local_index) {
    char file_out[100];
    int nintcf_g[nprocs], proc, i, **l2g_g, ncells;
    double **var_g, *var_cummulated;
    MPI_Status status;
    
    //perform communication to exchange data
    if (myrank == 0) {
        //TODO:analyze performance of this block, possibly create a single for loop
        //collect nintcf
        nintcf_g[0] = nintcf;
        for (proc=1; proc<nprocs; ++proc) {
            MPI_Recv(&nintcf_g[proc], 1, MPI_INT, proc, POSL_MPI_TAG_NINTCF, MPI_COMM_WORLD, &status);
        }
        
        //init global var array
        if ((var_g = (double**) malloc(nprocs*sizeof(double*))) == NULL) {
            fprintf(stderr, "malloc failed to allocate first dimension of var_g");
        }
        for (proc=0; proc<nprocs; ++proc) {
            if ((var_g[proc] = (double *) malloc((nintcf_g[proc]+1)*sizeof(double))) == NULL) {
                fprintf(stderr, "malloc failed to allocate second dimension of var_g\n");
            }
        }
        
        //collect var
        memcpy(var_g[0], var, (nintcf+1)*sizeof(double));
        for (proc=1; proc<nprocs; ++proc) {
            MPI_Recv(var_g[proc], (nintcf_g[proc]+1), MPI_DOUBLE, proc, POSL_MPI_TAG_VAR, 
                    MPI_COMM_WORLD, &status);
        }
        
        //init l2g_g
        if ((l2g_g = (int**) malloc(nprocs*sizeof(int*))) == NULL) {
            fprintf(stderr, "malloc failed to allocate first dimension of l2g_g");
        }
        for (proc=0; proc<nprocs; ++proc) {
            if ((l2g_g[proc] = (int *) malloc((nintcf_g[proc]+1)*sizeof(int))) == NULL) {
                fprintf(stderr, "malloc failed to allocate second dimension of l2g_g\n");
            }
        }
        
        //collect l2g_g
        memcpy(l2g_g[0], local_global_index, (nintcf+1)*sizeof(int));
        for (proc=1; proc<nprocs; ++proc) {
            MPI_Recv(l2g_g[proc], (nintcf_g[proc]+1), MPI_INT, proc, POSL_MPI_TAG_L2G, 
                    MPI_COMM_WORLD, &status);
        }
        
        //count total number of cells
        ncells = 0;
        for (proc=0; proc<nprocs; ++proc) {
            ncells += nintcf_g[proc]+1;
        }
        
        //build the final var
        if ((var_cummulated = (double *) malloc(ncells*sizeof(double))) == NULL) {
                fprintf(stderr, "malloc failed to allocate second dimension of var_cummulated\n");
        }
        for (proc=0; proc<nprocs; ++proc) {
            for (i=0; i<(nintcf_g[proc]+1); ++i){
                var_cummulated[l2g_g[proc][i]] = var_g[proc][i];
            }
        }
    } else {
        //send nintcf
        MPI_Send(&nintcf, 1, MPI_INT, 0, POSL_MPI_TAG_NINTCF, MPI_COMM_WORLD);
        
        //send var
        MPI_Send(var, nintcf+1, MPI_DOUBLE, 0, POSL_MPI_TAG_VAR, MPI_COMM_WORLD);
        
        //send l2g
        MPI_Send(local_global_index, nintcf+1, MPI_INT, 0, POSL_MPI_TAG_L2G, MPI_COMM_WORLD);
    }

    //perform stats printout
    if (myrank == 0) {
        sprintf(file_out, "%s_summary.out", file_in);
        int status = store_simulation_stats(file_in, file_out, nintci, ncells-1, var_cummulated, 
                total_iters, residual_ratio);
        if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
        
        
        //free memory
        for (proc=0; proc<nprocs; ++proc) {
            free(var_g[proc]);
            free(l2g_g[proc]);
        }
        free(var_g);
        free(l2g_g);
        free(var_cummulated);
    }
}
