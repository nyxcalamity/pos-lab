/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL 
 * computational arrays
 *
 * @author V. Petkov, A. Berariu, Denys Korzh, Denys Sobchyshak
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <papi.h>
#include "metis.h"

#include "initialization.h"
#include "initialization_algorithms.h"
#include "util_read_files.h"
#include "util_write_files.h"
#include "test_functions.h"
#include "util_processors.h"
#include "posl_definitions.h"
#include "util_errors.h"


int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank, 
        int* nintci, int* nintcf, int* nextci, int* nextcf, int*** lcc, double** bs, double** be,
        double** bn, double** bw, double** bl, double** bh, double** bp, double** su, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc, double** cnorm, 
        int** local_global_index, int** global_local_index, int *nghb_cnt, int** nghb_to_rank, 
        int** send_cnt, int*** send_lst,  int **recv_cnt, int*** recv_lst) {
    int input_key, part_key, read_key;
    process_cl(file_in, part_type, read_type, &input_key, &part_key, &read_key);
    /********** START INITIALIZATION **********/
    //TODO:optimize inits and var names
    // Used by metis function(gives us information to which process belongs our cell)
    int *partitioning_map, i=0;
    /** Simulation global variables which are read by first process to pass needed data
     * to other processes or which is needed by METIS */
    int nintci_g, nintcf_g;    /// internal cells start and end index
    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci_g, nextcf_g;
    int **lcc_g;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs_g, *be_g, *bn_g, *bw_g, *bh_g, *bl_g;
    double *bp_g;    /// Pole coefficient
    double *su_g;    /// Source values
    /** Geometry data */
    int points_count_g;    /// total number of points that define the geometry
    int** points_g;    /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems_g;    /// definition of the cells using their nodes (points) - each cell has 8 points
    /** Partitioning data */
    /// Metis will compute amount of cells for each process(classical example we will set it "manually")
    /// we can not be sure that all processors have the same amount of cells
    int int_cells_per_proc[nprocs];
    /// we can not be sure that all processors have the same amount of cells
    int extcell_per_proc[nprocs];
    int** local_global_index_g;

    int f_status = read_init_data(file_in, read_key, myrank,  &nintci_g, &nintcf_g, &nextci_g, 
            &nextcf_g, &lcc_g, &bs_g, &be_g, &bn_g, &bw_g, &bl_g, &bh_g, &bp_g, &su_g, &points_count_g, 
            &points_g, &elems_g);
    check_status(myrank, f_status, "Reading initial data failed");

    f_status = partition(part_key, read_key, myrank, nprocs, nintci_g, nintcf_g, nextci_g, nextcf_g, 
            &*nintci, &*nintcf, &*nextci, &*nextcf, lcc_g, points_count_g, points_g, elems_g, 
            int_cells_per_proc, extcell_per_proc, local_global_index_g, &*local_global_index, 
            &partitioning_map);
    check_status(myrank, f_status, "Domain partitioning failed");
    
    bcast_partitioning(read_key, myrank, &partitioning_map, &nintci_g, &nintcf_g, &nextci_g, &nextcf_g);
    
    f_status = allocate_lcc_elems_points(read_key, myrank, nprocs, nintci, nintcf, nextci, &*lcc, 
            &*points_count, &*points, &*elems, &*local_global_index, points_count_g, int_cells_per_proc);
    check_status(myrank, f_status, "Allocating elements failed");
    
    f_status = fill_l2g(read_key, myrank, nprocs, *nintcf, &*local_global_index, &local_global_index_g, 
            partitioning_map, nintcf_g-nintci_g+1, int_cells_per_proc);
    check_status(myrank, f_status, "Filling l2g failed");
    
    f_status = fill_lcc_elems_points(read_key, myrank, nprocs, *nintci, *nintcf, *lcc, *points_count, 
            *points, *elems, *local_global_index, local_global_index_g, lcc_g, points_count_g, points_g, 
            &elems_g, int_cells_per_proc);
    check_status(myrank, f_status, "Filling lcc and etc failed");
    
    f_status = build_lists_g2l_next(nprocs, myrank, partitioning_map, nintcf_g, nextcf_g, &*nintcf, &*nextcf, 
            &*lcc, &*local_global_index, &*global_local_index, &*nghb_cnt, &*nghb_to_rank, 
            &*recv_cnt, &*recv_lst);
    check_status(myrank, f_status, "Building g2l lists failed");
    
    f_status = allocate_send_lists(myrank, &*nghb_cnt, &*nghb_to_rank, &*send_cnt, &*send_lst, &*recv_cnt);
    check_status(myrank, f_status, "Allocating send lists failed");
    
    exchange_lists(myrank, &*nghb_cnt, &*nghb_to_rank, &*send_cnt, &*send_lst, &*recv_cnt, &*recv_lst);
    
    f_status = allocate_boundary_coef(nextcf, &*bs,&*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su);
    check_status(myrank, f_status, "Allocating boundary coefficients failed");
    
    f_status = fill_boundary_coef(read_key, myrank, nprocs, *nintci, *nintcf, *nextci, *nextcf, *bs, 
            *be, *bn, *bw, *bl, *bh, *bp, *su, *local_global_index, local_global_index_g, 
            &bs_g, &be_g, &bn_g, &bw_g, &bl_g, &bh_g, &bp_g, &su_g, int_cells_per_proc);
    check_status(myrank, f_status, "Filling boundary coefficients failed");

    // Check LCC
    if(OUTPUT_LCC_G) {
        if(myrank==0) {
            for(i=0;i<nintcf_g+1;++i) {
                printf("i%-6d,p%d, %-10d %-10d  %-10d %-10d %-10d %-10d\n",
                        i,partitioning_map[i],lcc_g[i][0],lcc_g[i][1],lcc_g[i][2],lcc_g[i][3],lcc_g[i][4],
                        lcc_g[i][5]);
            }
        }
    }
    if (OUTPUT_NINTCF_NINTCE) {
        printf("rank%d,*nintci=%d,*nintcf=%d, *nextci=%d, *nextcf=%d\n", myrank,*nintci,*nintcf,*nextci,*nextcf);
    }
    // End check lcc
    
    if (OUTPUT_LCC) {
        if (myrank==0) {
            for (i=0;i<(*nintcf)+1;++i) {
                printf("i%-6d, %d\n", i, (*global_local_index)[i]);
            }
        }
    }

    *var = (double*) calloc(sizeof(double), (*nextcf+1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf+1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf+1));

    // initialize the arrays
    for (i=0; i<=10; i++) {
        (*cnorm)[i] = 1.0;
    }

    for (i=(*nintci); i<=(*nintcf); i++) {
        (*var)[i] = 0.0;
    }

    for (i=(*nintci); i<=(*nintcf); i++) {
        (*cgup)[i] = 1.0/((*bp)[i]);
    }

    for (i=(*nextci); i<=(*nextcf); i++) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bh)[i] = 0.0;
        (*bl)[i] = 0.0;
    }

    // VTK check
    if (OUTPUT_VTK) {
        vtk_check_lists(file_in, myrank,
                *local_global_index, (*nintcf-*nintci+1),
                *nghb_cnt, *nghb_to_rank, *send_cnt, *send_lst,
                *recv_cnt, *recv_lst, OUTPUT_VTK);
        vtk_check_neighbour(file_in, myrank,
                        *local_global_index, (*nintcf-*nintci+1),
                        *nghb_cnt, *nghb_to_rank, *send_cnt, *send_lst,
                        *recv_cnt, *recv_lst, OUTPUT_VTK, VTK_NEIGHBOUR);
    }
    
    //convert global indexes
    converte_global2local_idx(myrank, *global_local_index, *nintci, *nintcf, *lcc, *nghb_cnt,
            *send_cnt, *send_lst, *recv_cnt, *recv_lst);

    // Free
    //TODO:code proper memory freeing
    if((read_key == POSL_INIT_ONE_READ && myrank == 0) || read_key == POSL_INIT_ALL_READ) {
//        free(local_global_index_g);
        for (i = 0; i < (*nintci + 1); i++) {
            free(lcc_g[i]);
        }
        free(lcc_g);
        free(su_g);
        free(bp_g);
        free(bh_g);
        free(bl_g);
        free(bw_g);
        free(bn_g);
        free(be_g);
        free(bs_g);
        free(elems_g);
        free(*global_local_index);
    }
    
    if (DEBUG_ENABLED) {
        log_dbg("Initialization phase complete on process #%d", myrank);
    }
    /********** END INITIALIZATION **********/
    for (i=0; i<*nghb_cnt; ++i) {
        write_pstats_communication(input_key, part_key, myrank, nprocs, *nghb_cnt, i, 
                *send_cnt, *send_lst, *recv_cnt, *recv_lst);
    }
    return 0;
}
