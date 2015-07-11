#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "initialization.h"
#include "compute_solution.h"
//start_of_student_code-------------------------------------------------------------------------------------------------
#include "final_output_vtk.h"
#include <papi.h>
//end_of_student_code---------------------------------------------------------------------------------------------------
#include "finalization.h"

int main(int argc, char *argv[]) {
//start_of_student_code-------------------------------------------------------------------------------------------------
    if (argc != 4) {
        fprintf(stderr, "Wrong number of arguments:gccg <format> <input file> <output prefix>\n");
        abort();
    }
//end_of_student_code---------------------------------------------------------------------------------------------------
    int i;

    const int max_iters = 10000;    /// maximum number of iteration to perform

    /** Simulation parameters parsed from the input datasets */
    int nintci, nintcf;    /// internal cells start and end index
    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci, nextcf;
    int **lcc;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bh, *bl;
    double *bp;    /// Pole coefficient
    double *su;    /// Source values

    double residual_ratio;    /// the ratio between the reference and the current residual
    double *var;    /// the variation vector -> keeps the result in the end

    /** Additional vectors required for the computation */
    double *cgup, *oc, *cnorm;

//start_of_student_code-------------------------------------------------------------------------------------------------
    char *file_out_prefix = argv[3];
    char *file_in = argv[2];
    char *file_format = argv[1];
    if (strcmp(file_format, "bin") != 0 && strcmp(file_format, "text") != 0) {
        fprintf(stderr, "Wrong file format(should be bin or text)\n");
        abort();
    }
//end_of_student_code---------------------------------------------------------------------------------------------------
    /********** START INITIALIZATION **********/
    // read-in the input file
    int init_status = initialization(file_format, file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
            &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &var, &cgup, &oc,
            &cnorm);

    if (init_status != 0) {
        fprintf(stderr, "Failed to initialize data!\n");
        abort();
    }


    /********** END INITIALIZATION **********/
//start_of_student_code-------------------------------------------------------------------------------------------------
    const int EVENT_NUM=4;
    //int events[] = {PAPI_L2_TCM, PAPI_L2_TCA, PAPI_L3_TCM, PAPI_L3_TCA};
    int event_set = PAPI_NULL;
    long_long start_usec, end_usec, values[EVENT_NUM];

    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) exit(1);
    if (PAPI_create_eventset(&event_set) != PAPI_OK){
        printf("Failed to create eventset");
	exit(1);
    }
//    if (PAPI_add_events(event_set, events, EVENT_NUM) != PAPI_OK){
//         printf("Failed to add events");
//	   exit(1);
//    }
    if (PAPI_add_event(event_set, PAPI_DP_OPS) != PAPI_OK){
         printf("Failed to add events");
	exit(1);
    }

    if (PAPI_start(event_set) != PAPI_OK) {
        printf("Failed to start eventset");
	exit(1);
    }
    start_usec = PAPI_get_virt_usec();
//end_of_student_code---------------------------------------------------------------------------------------------------
    /********** START COMPUTATIONAL LOOP **********/
    int total_iters = compute_solution(max_iters, nintci, nintcf, nextcf, lcc, bp, bs, bw, bl, bn,
                                       be, bh, cnorm, var, su, cgup, &residual_ratio);
    /********** END COMPUTATIONAL LOOP **********/
//start_of_student_code-------------------------------------------------------------------------------------------------
    end_usec = PAPI_get_virt_usec();
    if(PAPI_read(event_set, values) != PAPI_OK){
        printf("Failed to read values");
	exit(1);
    }
    if(PAPI_stop(event_set, values) != PAPI_OK){
        printf("Failed to stop counters");
	exit(1);
    }

    printf("Computation time (secs):\t%f\n", (double)(end_usec-start_usec)/(double)1000000);
    
//    printf("L2 cache miss rate: \t%.3f%%\nL3 cache miss rate: \t%.3f%%\n",
//            (double)values[0]/(double)values[1], (double)values[2]/(double)values[3]);
    printf("MFLOPS:\t%f\n", values[0]/(double)(end_usec-start_usec));
//end_of_student_code---------------------------------------------------------------------------------------------------
    /********** START FINALIZATION **********/
    finalization(file_in, total_iters, residual_ratio, nintci, nintcf, var, cgup, su);
//start_of_student_code-------------------------------------------------------------------------------------------------
    finalOutputVTK(file_out_prefix, nintci, nintcf, lcc, var, cgup, su);
//end_of_student_code---------------------------------------------------------------------------------------------------
    /********** END FINALIZATION **********/
    free(cnorm);
    free(var);
    free(cgup);
    free(su);
    free(bp);
    free(bh);
    free(bl);
    free(bw);
    free(bn);
    free(be);
    free(bs);

for ( i = 0; i < 6; ++i)
{
   printf("lcc[1][%d] = %d\n", i, lcc[1][i] );
}

   for ( i = nintci; i <= nintcf; i++ ) {
        free(lcc[i]);
    }
   free(lcc);
 
    return 0;
}

