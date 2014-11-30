/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 13-Nov-2014
 * @author V. Petkov, A. Berariu
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

int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank, 
        int* nintci, int* nintcf, int* nextci, int* nextcf, int*** lcc, double** bs, double** be,
        double** bn, double** bw, double** bl, double** bh, double** bp, double** su, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc, double** cnorm, 
        int** local_global_index, int** global_local_index, int *nghb_cnt, int** nghb_to_rank, 
        int** send_cnt, int*** send_lst,  int **recv_cnt, int*** recv_lst) {
    long_long start_usec, end_usec;
    int input_key, part_key, read_key;
    process_cl(file_in, part_type, read_type, &input_key, &part_key, &read_key);
    start_usec = PAPI_get_virt_usec();
    
    /********** START INITIALIZATION **********/
    int m,n;
    //FIXME:optimize inits and var names
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
    //TODO:externalize error checking
    if (f_status != 0) {
        return f_status;
    }

    f_status = partition(part_key, read_key, myrank, nprocs, nintci_g, nintcf_g, nextci_g, nextcf_g, 
            &*nintci, &*nintcf, &*nextci, &*nextcf, lcc_g, points_count_g, points_g, elems_g, 
            int_cells_per_proc, extcell_per_proc, local_global_index_g, &*local_global_index, 
            &partitioning_map);
    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
    }
    
    if (DEBUG_OUTPUT_PARTITIONING && read_key == POSL_INIT_ONE_READ && myrank == 0) {
        printf("Int_cells_per_proc at proc #%d\n", myrank);
        for (m=0; m<nprocs; ++m) {
            printf("%d ", int_cells_per_proc[m]);
        }
        printf("\n");
    }
    bcast_partitioning(read_key, myrank, &partitioning_map, &nintci_g, &nintcf_g, &nextci_g, &nextcf_g);
    
    if (DEBUG_OUTPUT_PARTITIONING && read_key == POSL_INIT_ONE_READ) {
        printf("partition_map at proc #%d\n", myrank);
        for (m=0; m<30; ++m) {
            printf("%d ", partitioning_map[m]);
        }
        printf("\n");
        for (m=nextci_g-1; m>(nextci_g-30); --m) {
            printf("%d ", partitioning_map[m]);
        }
        printf("\n");
    }
    
    
    f_status = allocate_lcc_elems_points(read_key, myrank, nprocs, nintci, nintcf, &*lcc, &*points_count, 
            &*points, &*elems, &*local_global_index, points_count_g, int_cells_per_proc);
    
//    printf("#%d number of global elements = %d\n", myrank, nintcf_g-nintci_g+1);
    fill_l2g(read_key, myrank, nprocs, *nintcf, &*local_global_index, &local_global_index_g, 
            partitioning_map, nintcf_g-nintci_g+1, int_cells_per_proc);
    
    if (DEBUG_OUTPUT_L2G_G) {
        printf("local2global index at proc #%d\n", myrank);
        for (m=0; m<20; ++m) {
            printf("%d ", local_global_index[m]);
        }
        printf("\n");
    }
    
    f_status = fill_lcc_elems_points(read_key, myrank, nprocs, *nintci, *nintcf, *lcc, *points_count, 
            *points, *elems, *local_global_index, local_global_index_g, lcc_g, points_count_g, points_g, &elems_g, 
            int_cells_per_proc);
    printf("Rank #%d got past lcc filling\n", myrank);
    
    build_lists_g2l_next(nprocs, myrank, partitioning_map, nintcf_g, nextcf_g, &*nintcf, &*nextcf, 
            &*lcc, &*local_global_index, &*global_local_index, &*nghb_cnt, &*nghb_to_rank, 
            &*recv_cnt, &*recv_lst);
    printf("Rank #%d got past list building\n", myrank);
    
    allocate_send_lists(myrank, &*nghb_cnt, &*nghb_to_rank, &*send_cnt, &*send_lst, &*recv_cnt);
    printf("Rank #%d got past send list building\n", myrank);
    
    exchange_lists(myrank, &*nghb_cnt, &*nghb_to_rank, &*send_cnt, &*send_lst, &*recv_cnt, &*recv_lst);
    printf("Rank #%d got past list exchange\n", myrank);
    
    f_status = allocate_boundary_coef(nextcf, &*bs,&*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su);
    printf("Rank #%d got past boundary coef allocation\n", myrank);
    
    f_status = fill_boundary_coef(read_key, myrank, nprocs, *nintci, *nintcf, *nextci, *nextcf, *bs, 
            *be, *bn, *bw, *bl, *bh, *bp, *su, *local_global_index, local_global_index_g, 
            &bs_g, &be_g, &bn_g, &bw_g, &bl_g, &bh_g, &bp_g, &su_g, int_cells_per_proc);
    printf("Rank #%d got past boundary coef filling\n", myrank);

    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
    }

    // Check LCC
    if(OUTPUT_LCC_G) {
        if(myrank==0) {
            for(i=0;i<nintcf_g+1;++i) {
                printf("i%-6d,p%d, %-10d %-10d  %-10d %-10d %-10d %-10d\n",
                        i,partitioning_map[i],lcc_g[i][0],lcc_g[i][1],lcc_g[i][2],lcc_g[i][3],lcc_g[i][4],
                        lcc_g[i][5]);
//                printf("i%-6d, %d\n", i, local_global_index_g[i]);
            }
        }
    }
    if (OUTPUT_NINTCF_NINTCE) {
        printf("rank%d,*nintcf=%d, *nextcf=%d\n", myrank,*nintcf,*nextcf);
    } // End check lcc
    
    if (OUTPUT_LCC) {
        if (myrank==0) {
            for (i=0;i<(*nintcf)+1;++i) {
//                printf("i%-6d, %-10d %-10d  %-10d %-10d %-10d %-10d\n",
//                        i,(*lcc)[i][0],(*lcc)[i][1],(*lcc)[i][2],(*lcc)[i][3],(*lcc)[i][4],(*lcc)[i][5]);
//                printf("i%-6d, %d\n", i, (*local_global_index)[i]);
                printf("i%-6d, %d\n", i, (*global_local_index)[i]);
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
    if (OUTPUT_VTK) {
        f_status = vtk_check(file_in, myrank, *nintci, *nintcf, *su, *cgup, *points_count, *points,
                *elems, *local_global_index, (*nintcf-*nintci+1));
        char szFileName[80];
        sprintf(szFileName, "%s%s.receive.rank%i.vtk", "out/", "test", myrank);
        test_distribution(file_in, szFileName, (*recv_lst)[0], (*recv_cnt)[0], *cgup);
        sprintf(szFileName, "%s%s.send.rank%i.vtk", "out/", "test", myrank);
        test_distribution(file_in, szFileName, (*send_lst)[0], (*send_cnt)[0], *cgup);
    }
    //TODO:externalize error checking
    if (f_status != 0){
        return f_status;
    }

    // Free
    //FIXME:code proper memory freeing
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
        // It is not freed in gccg.c
        free(*global_local_index);
    }
//    printf("[INFO] Completed initialization on task #%d\n", myrank);
    /********** END INITIALIZATION **********/
    
    end_usec = PAPI_get_virt_usec();
    write_pstats_exectime(input_key, part_key, read_key, myrank, (double)(end_usec-start_usec));
    write_pstats_partition(input_key, part_key, myrank, int_cells_per_proc[myrank], extcell_per_proc[myrank]);

    printf("Rank #%d completed execution\n", myrank);
    return 0;
}
