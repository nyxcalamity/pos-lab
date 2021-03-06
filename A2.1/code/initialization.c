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


int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index) {
    long_long start_usec, end_usec;
    int input_key, part_key, read_key;
    process_cl(file_in, part_type, read_type, &input_key, &part_key, &read_key);
//    start_usec = PAPI_get_virt_usec();
    /********** START INITIALIZATION **********/
    //FIXME:optimize inits and var names
    // Used by metis function(gives us information to which process belongs our cell)
    int *metis_idx, i=0;
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
    int intcell_per_proc[nprocs];
    /// we can not be sure that all processors have the same amount of cells
    int extcell_per_proc[nprocs];
    int* local_global_index_g;

    // Initialize and read all global variables or read data for geometry
    int f_status = read_global_data_or_geometry(file_in, read_type, myrank,  &nintci_g, &nintcf_g, 
            &nextci_g, &nextcf_g, &lcc_g, &bs_g, &be_g, &bn_g, &bw_g, &bl_g, &bh_g, &bp_g, &su_g, 
            &points_count_g, &points_g, &elems_g);
    //TODO:externalize error checking
    if (f_status != 0) {
        return f_status;
    }
    // METIS/NUMBER OF CELLS FOR EACH PROCESSOR/LOCAL TO GLOBAL
    f_status = compute_metis(part_type,read_type, myrank, nprocs, nintci_g, nintcf_g, nextci_g, 
            nextcf_g, &*nintci, &*nintcf, &*nextci, &*nextcf, lcc_g, points_count_g, points_g, elems_g,
            intcell_per_proc, extcell_per_proc, &local_global_index_g, &metis_idx);
    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
    }
    if(OUTPUT_LCC_G) {
        if(myrank==0) {
            for(i=0;i<nintcf_g+1;++i) {
                printf("i%-6d,p%d, %-10d %-10d  %-10d %-10d %-10d %-10d\n",
                        i,metis_idx[i],lcc_g[i][0],lcc_g[i][1],lcc_g[i][2],lcc_g[i][3],lcc_g[i][4],lcc_g[i][5]);
//                printf("i%-6d, %d\n", i, local_global_index_g[i]);
            }
        }
    }
    if (OUTPUT_NINTCF_NINTCE) {
        printf("rank%d,intcell_per_proc[end]=%d, extcell_per_proc[end]=%d\n",
                myrank,intcell_per_proc[nprocs-1],extcell_per_proc[nprocs-1]);
        printf("rank%d,intcell_per_proc[0]=%d, extcell_per_proc[0]=%d\n",
                myrank,intcell_per_proc[0],extcell_per_proc[0]);
    }
    // ALLOCATE lcc, bs, be, bn, bw, bl, bh, bp, su, points, local_global_index
    f_status = allocate_local_variables(read_type, myrank, nprocs, nintci, nintcf, nextci, nextcf,
            &*lcc, &*bs,&*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count, &*points, &*elems, 
            &*local_global_index, intcell_per_proc, extcell_per_proc, local_global_index_g, points_count_g);
    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
    }
    f_status = send_or_read_data(read_type, myrank, nprocs, *nintci, *nintcf, *nextci, *nextcf, *lcc, 
            *bs, *be, *bn, *bw, *bl, *bh, *bp, *su, *points_count, *points, *elems, 
            *local_global_index, intcell_per_proc, extcell_per_proc, nintci_g, nintcf_g, nextci_g, 
            nextcf_g, lcc_g, &bs_g, &be_g, &bn_g, &bw_g, &bl_g, &bh_g, &bp_g, &su_g, points_count_g, 
            points_g, &elems_g, local_global_index_g);

    if (OUTPUT_LCC) {
        if (myrank==1) {
            for(i=0;i<(*nintcf)+1;++i) {
//                printf("i%-6d, %-10d %-10d  %-10d %-10d %-10d %-10d\n",
//                        i,(*lcc)[i][0],(*lcc)[i][1],(*lcc)[i][2],(*lcc)[i][3],(*lcc)[i][4],(*lcc)[i][5]);
                printf("i%-6d, %d\n", i, (*local_global_index)[i]);
            }
        }
    }

    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
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
    f_status = vtk_check(file_in, myrank, *nintci, *nintcf, *su, *cgup, *points_count, *points,
            *elems, *local_global_index, (*nintcf-*nintci+1));
    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
    }

    // Free
    //FIXME:code proper memory freeing
    if((read_key == POSL_INIT_ONE_READ && myrank == 0) || read_key == POSL_INIT_ALL_READ) {
        free(local_global_index_g);
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
    }
    printf("[INFO] Completed initialization on task #%d\n", myrank);
    /********** END INITIALIZATION **********/
//    end_usec = PAPI_get_virt_usec();
//    write_pstats_exectime(input_key, part_key, read_key, myrank, (double)(end_usec-start_usec));
//    write_pstats_partition(input_key, part_key, myrank, intcell_per_proc[myrank], extcell_per_proc[myrank]);
    return 0;
}
