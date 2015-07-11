#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#include "initialization.h"
#include "compute_solution.h"
#include "finalization.h"
#include "test_functions.h"
#include "util_errors.h"
#include "util_processors.h"


int main(int argc, char *argv[]) {
    int my_rank, num_procs, i;

    const int max_iters = 600;    /// maximum number of iteration to perform

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

    /** Geometry data */
    int points_count;    /// total number of points that define the geometry
    int** points;    /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems;    /// definition of the cells using their nodes (points) - each cell has 8 points

    /** Mapping between local and remote cell indices */
    int* local_global_index;    /// local to global index mapping
    int* global_local_index;    /// global to local index mapping
  

    /** Lists for neighbouring information */
    int nghb_cnt = 0;    /// total number of neighbors of the current process
    int *nghb_to_rank;  /// mapping of the neighbour index to the corresponding process rank
    int *send_cnt;    /// number of cells to be sent to each neighbour (size: nghb_cnt)
    int **send_lst;    /// lists of cells to be sent to each neighbour (size: nghb_cnt x send_cnt[*])
    int *recv_cnt;    /// number of cells to be received from each neighbour (size: nghb_cnt)
    int **recv_lst;    /// lists of cells to be received from each neighbour (size: nghb_cnt x recv_cnt[*])
    
    // time measurement variables
    double start_time, end_time, min_start_time, max_end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //process command line arguments
    if (argc < 4) {
        log_err("Usage: ./gccg <input_file> <partition_type> <algorithm_type>");
        MPI_Abort(MPI_COMM_WORLD, POSL_ERROR);
    }

    char *file_in = argv[1];
    char *part_type = argv[2];
    char *read_type = argv[3];
    int input_key, part_key, read_key;
    process_cl(file_in, part_type, read_type, &input_key, &part_key, &read_key);
    
    if (part_key == POSL_ERROR) {
        log_err("Wrong partition type selected. Valid values are classic, nodal and dual");
        MPI_Abort(MPI_COMM_WORLD, POSL_ERROR);
    }
    
    if (read_key == POSL_ERROR) {
        log_err("Wrong read-in algorithm selected. Valid values are oneread and allread");
        MPI_Abort(MPI_COMM_WORLD, POSL_ERROR);
    }

    /********** START INITIALIZATION **********/
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    int init_status = initialization(file_in, input_key, part_key, read_key, num_procs, my_rank,
                                     &nintci, &nintcf, &nextci, &nextcf, 
                                     &lcc, &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, 
                                     &points_count, &points, &elems, &var, &cgup, &oc, &cnorm, 
                                     &local_global_index, &global_local_index,
                                     &nghb_cnt, &nghb_to_rank, 
                                     &send_cnt, &send_lst, &recv_cnt, &recv_lst);
    end_time = MPI_Wtime();
    //find min and max times
    MPI_Reduce(&start_time, &min_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //print status message
    if (my_rank == 0) {
        log_info("Initialization ET (secs): %f", max_end_time-min_start_time);
    } 
    /** LOCAL DATA FROM HERE ON **/
    // at this point, all initialized vectors should contain only the locally needed data
    // and all variables representing the number of elements, cells, points, etc. should 
    // reflect the local setup, e.g. nintcf-nintci+1 is the local number of internal cells
    if (init_status != 0) {
        log_err("Failed to initialize data!");
        MPI_Abort(MPI_COMM_WORLD, POSL_ERROR);
    }
    /********** END INITIALIZATION **********/

    /********** START COMPUTATIONAL LOOP **********/
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    int total_iters = compute_solution(num_procs, my_rank, max_iters, nintci, nintcf, nextcf, 
                    lcc, bp, bs, bw, bl, bn, be, bh,
                     cnorm, var, su, cgup, &residual_ratio,
                     local_global_index, global_local_index, nghb_cnt, 
                     nghb_to_rank, send_cnt, send_lst, recv_cnt, recv_lst);
    end_time = MPI_Wtime();
    //find min and max times
    MPI_Reduce(&start_time, &min_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //print status message
    if (my_rank == 0) {
        log_info("Computation ET (secs): %f", max_end_time-min_start_time);
    }
    /********** END COMPUTATIONAL LOOP **********/
    
    /********** START FINALIZATION **********/
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();
    finalization(file_in, num_procs, my_rank, total_iters, residual_ratio, nintci, nintcf, var, 
            local_global_index, global_local_index);
    end_time = MPI_Wtime();
    //find min and max times
    MPI_Reduce(&start_time, &min_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &max_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //print status message
    if (my_rank == 0) {
        log_info("Finalization ET (secs): %f", max_end_time-min_start_time);
    }
    /********** END FINALIZATION **********/

    // cleanup allocated memory
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
    free(elems);

    for ( i = 0; i < nintcf + 1; i++ ) {
        free(lcc[i]);
    }
    free(lcc);

    for ( i = 0; i < points_count; i++ ) {
        free(points[i]);
    }
    free(points);

    free(nghb_to_rank);
    
    free(send_cnt);
    for ( i = 0; i < nghb_cnt; i++ ){
        free(send_lst[i]);
    }
    free(send_lst);

    free(recv_cnt);
    for ( i = 0; i < nghb_cnt; i++ ){
        free(recv_lst[i]);
    }
    free(recv_lst);

    MPI_Finalize();
    
    return 0;
}
