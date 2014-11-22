#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "metis.h"

#include "initialization_algorithms.h"
#include "util_read_files.h"
#include"posl_definitions.h"


int read_global_data_or_geometry(char* file_in,char* read_type, int myrank, int *nintci, int *nintcf, 
        int *nextci, int *nextcf, int ***lcc, double **bs, double **be, double **bn, double **bw, 
        double **bl, double **bh, double **bp, double **su, int* points_count, int***points, 
        int** elems) {
    int f_status=0;
    //TODO:replace strcmp with definitions
    if (!strcmp(read_type, "oneread")) {
        if (myrank == 0) {
            // read-in the input file
            f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                    &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count, &*points, &*elems);
        }
    } else {
        //TODO:implemnete sheer geometry reading
        f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs, 
                &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
                &*points, &*elems);
    }
    if (f_status != 0) {
        return f_status;
    } else {
        return 0;
    }
}


//FIXME:rename this function
int compute_metis(char* part_type, char* read_type, int myrank, int nprocs, int nintci_g, 
        int nintcf_g, int nextci_g, int nextcf_g, int *nintci, int *nintcf, int *nextci, int *nextcf, 
        int **lcc_g, int points_count_g, int**points_g, int* elems_g, int *intcell_per_proc, 
        int *extcell_per_proc, int** local_global_index_g, int** local_global_index, int **metis_idx) {
    int i=0;
    idx_t nelems, nnodes, ncommon, nparts, objval;
    idx_t *elem_ptr, *elem_idx, *elem_part, *node_part;
    nelems = nintcf_g-nintci_g+1;
    *metis_idx = (int *) calloc(sizeof(int), (nintcf_g-nintci_g+1));
    //TODO: replace strcmp with definitions
    if ((!strcmp(read_type, "oneread") && (myrank == 0))) {

    } else if (!strcmp(read_type, "allread")) {
        *nintci = 0;
        *nintcf = 0;
        if (!strcmp(part_type, "classic")) {
            if(myrank == nprocs-1) {
                *nextci = nelems-(nprocs-1)*((nelems+(nprocs-1))/nprocs);
                *nintcf = *nextci - 1;
            } else {
                *nextci = (nelems+(nprocs-1))/nprocs;
                *nintcf = *nextci - 1;
            }
            for (i=0; i<nelems; ++i) {
                (*metis_idx)[i] = i / ((nelems+(nprocs-1))/nprocs);
            }
        } else {
            nnodes = points_count_g;
            ncommon = 4;
            nparts = nprocs;

            elem_ptr = (idx_t *) calloc(nelems+1, sizeof(idx_t));
            elem_idx = (idx_t *) calloc(nelems*8, sizeof(idx_t));
            elem_part = (idx_t *) calloc(nelems, sizeof(idx_t));
            node_part = (idx_t *) calloc(nnodes, sizeof(idx_t));

            for (i=0; i<(nelems+1); i++) {
                elem_ptr[i] = 8*i;
            }
            for (i=0; i<(nelems*8); i++) {
                elem_idx[i] = elems_g[i];
            }
            //TODO: replace strcmp with definitions
            if (!strcmp(part_type, "dual")) {
                METIS_PartMeshDual(&nelems, &nnodes, elem_ptr, elem_idx, NULL, NULL, &ncommon, &nparts, NULL, NULL,
                        &objval, elem_part, node_part);
            } else {
                METIS_PartMeshNodal(&nelems, &nnodes, elem_ptr, elem_idx, NULL, NULL, &nparts, NULL, NULL, &objval,
                        elem_part, node_part);
            }
            for (i=0; i<nelems; i++) {
                (*metis_idx)[i] = (int) elem_part[i];
            }
            // Calculate nintcf
            for (i=0; i<nelems; i++) {
                if (myrank == (*metis_idx)[i]) {
                    (*nintcf) += 1;
                }
            }
            // At the end we need to subtract one because nintcf is the last index not the amount of elements
            *nextci = (*nintcf)--;
        }
    }
    return 0;
}


int allocate_lcc_elems_points(char* read_type, int myrank, int nprocs, int *nintci, int *nintcf,
       int ***lcc, int* points_count, int*** points,
        int** elems, int **local_global_index, int points_count_g) {
    int i=0;
    MPI_Status status;
    //TODO: replace strcmp with definitions
    if (!strcmp(read_type, "oneread")) {

    } else {
        *points_count = points_count_g;
    }
    if ((*lcc = (int**) malloc(((*nintcf)+1)*sizeof(int*))) == NULL) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }
    for (i=0; i<(*nintcf)+1; i++) {
        if (((*lcc)[i] = (int *) malloc(6*sizeof(int))) == NULL) {
            fprintf(stderr, "malloc failed to allocate second dimension of lcc\n");
            return -1;
        }
    }
    if ((*elems = (int *) malloc(((*nintcf)+1)*8*sizeof(int))) == NULL) {
        fprintf(stderr, "malloc(elems) failed\n");
        return -1;
    }
    if ( (*points = (int **) calloc(*points_count, sizeof(int*))) == NULL) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }
    for ( i=0; i<*points_count; i++) {
        if (((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }
    return 0;
}


int fill_lcc_elems_points(char* read_type, int myrank, int nprocs, int nintci, int nintcf,
        int **lcc, int points_count, int** points, int* elems,
        int *local_global_index,  int **lcc_g,
        int points_count_g, int** points_g, int **elems_g) {
    int k=0, i=0;
    MPI_Status status;
    //TODO: replace strcmp with definitions
    if (!strcmp(read_type, "oneread")) {

    } else {
        for(i=nintci; i<nintcf+1; i++) {
            memcpy(lcc[i], lcc_g[local_global_index[i]], 6*sizeof(int));
            memcpy(&(elems[8*i]), &(*elems_g)[local_global_index[i]*8], 8*sizeof(int));
        }
        for (i=0; i<points_count; i++) {
            memcpy(points[i], points_g[i], 3*sizeof(int));
        }
    }
    return 0;
}


int allocate_boundary_coef(char* read_type, int myrank, int nprocs, int *nextcf, double **bs, double **be, double **bn, double **bw,
        double **bl, double **bh, double **bp, double **su) {
    if ((*bs = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*be = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bn = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bw = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bl = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bh = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bp = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*su = (double *) malloc(((*nextcf)+1)*sizeof(double))) == NULL) {
        printf("malloc() failed\n");
        return -1;
    }
    return 0;
}


int fill_boundary_coef(char* read_type, int myrank, int nprocs, int nintci, int nintcf, int nextci,
        int nextcf, double *bs, double *be, double *bn, double *bw, double *bl,
        double *bh, double *bp, double *su,
        int *local_global_index, double **bs_g, double **be_g,
        double **bn_g, double **bw_g, double **bl_g, double **bh_g, double **bp_g, double **su_g) {
    int k=0, i=0;
    MPI_Status status;
    //TODO: replace strcmp with definitions
    if (!strcmp(read_type, "oneread")) {

    } else {
        for(i=nintci; i<nintcf+1; i++) {
            bs[i] = (*bs_g)[local_global_index[i]];
            be[i] = (*be_g)[local_global_index[i]];
            bn[i] = (*bn_g)[local_global_index[i]];
            bw[i] = (*bw_g)[local_global_index[i]];
            bl[i] = (*bl_g)[local_global_index[i]];
            bh[i] = (*bh_g)[local_global_index[i]];
            bp[i] = (*bp_g)[local_global_index[i]];
            su[i] = (*su_g)[local_global_index[i]];
        }
    }
    return 0;
}


/**
 * Compute local_global_index in sorted order
 */
//FIXME: change name of this function + check on bad allocation
void fill_local_global_index(char* read_type, int myrank, int nintci, int nintcf,
        int** local_global_index, int *metis_idx, int nelems_g) {
    int i = 0, local_idx=0;
    if ((*local_global_index = (int *) malloc(((nintcf)+1)*sizeof(int))) == NULL) {
        fprintf(stderr, "malloc(local_global_index) failed\n");
//        return -1;
    }
    for (i=0; i<nelems_g; ++i) {
        if (metis_idx[i]==myrank) {
            (*local_global_index)[local_idx]=i;
            ++local_idx;
        }
    }
}


int sort_data_by_local_global_index(int nintci_g, int nintcf_g, int nextci_g, int nextcf_g, 
        double **bs_g, double **be_g, double **bn_g, double **bw_g, double **bl_g, 
        double **bh_g, double **bp_g, double **su_g, int **elems_g, int *local_global_index_g) {
    int i=0;
    // This two arrays are needed to send data in correct order
    int *tmp_elems;
    double *tmp_b;
    // Not to lose the pointer
    void *tmp;
    if ((tmp_elems = (int *) malloc((nintcf_g+1)*8*sizeof(int))) == NULL) {
        printf("malloc() of tmp_elems in send_or_read_data failed\n");
        return -1;
    }
    if ((tmp_b = (double *) malloc((nintcf_g+1)*sizeof(double))) == NULL) {
        printf("malloc() of tmp_b in send_or_read_data failed\n");
        return -1;
    }
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*bs_g)[local_global_index_g[i]];
        memcpy(&tmp_elems[8*i], &(*elems_g)[8*local_global_index_g[i]],8*sizeof(int));
    }
    //TODO:write a function to do the swap operations
    tmp = *elems_g;
    *elems_g = tmp_elems;
    tmp_elems = tmp;
    tmp=*bs_g;
    *bs_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*be_g)[local_global_index_g[i]];
    }
    tmp=*be_g;
    *be_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*bn_g)[local_global_index_g[i]];
    }
    tmp=*bn_g;
    *bn_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*bw_g)[local_global_index_g[i]];
    }
    tmp=*bw_g;
    *bw_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*bl_g)[local_global_index_g[i]];
    }
    tmp=*bl_g;
    *bl_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*bh_g)[local_global_index_g[i]];
    }
    tmp=*bh_g;
    *bh_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*bp_g)[local_global_index_g[i]];
    }
    tmp=*bp_g;
    *bp_g=tmp_b;
    tmp_b=tmp;
    for (i=0; i<nintcf_g+1; ++i) {
        tmp_b[i] = (*su_g)[local_global_index_g[i]];
    }
    tmp=*su_g;
    *su_g=tmp_b;
    tmp_b=tmp;
    free(tmp_b);
    free(tmp_elems);
    return 0;
}


void build_lists_g2l_next(char* part_type, char* read_type, int nprocs, int myrank,
        int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst,
        int **recv_cnt, int*** recv_lst) {
    if ((!strcmp(read_type, "oneread") && (myrank == 0)) || !strcmp(read_type, "allread")) {
        // TODO: do it!
        *nextcf = *nintcf + 10000;
    }
}
