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
        int *extcell_per_proc, int** local_global_index_g, int **metis_idx) {
    int start_idx=0, proc=0, i=0;
    idx_t nelems, nnodes, ncommon, nparts, objval;
    idx_t *elem_ptr, *elem_idx, *elem_part, *node_part;
    //TODO: replace strcmp with definitions
    if ((!strcmp(read_type, "oneread") && (myrank == 0)) || !strcmp(read_type, "allread")) {
        if ((*local_global_index_g = (int *) malloc((nintcf_g + 1) * sizeof(int))) == NULL) {
            fprintf(stderr, "malloc(local_global_index) failed\n");
            return -1;
        }
        *metis_idx = (int *) calloc(sizeof(int), (nintcf_g-nintci_g+1));
        //TODO: replace strcmp with definitions
        if (!strcmp(part_type, "classic")) {
            for (i=0; i<nprocs-1; i++) {
                //adjusted devision to account for leftover cells
                intcell_per_proc[i] = (nintcf_g+1+(nprocs-1))/nprocs;
            }
            //if our domain can't be divided in equal parts
            intcell_per_proc[nprocs-1] = nintcf_g+1-(nprocs-1)*((nintcf_g+1+(nprocs-1))/nprocs);
            for (i=0; i<(nintcf_g+1); i++) {
                (*local_global_index_g)[i] = i;
            }
            for (proc=0; proc<nprocs; ++proc) {
                for (i=0; i<intcell_per_proc[proc]; ++i) {
                    (*metis_idx)[start_idx+i] = proc;
                }
                start_idx += intcell_per_proc[proc];
            }
        } else {
            nelems = nintcf_g-nintci_g+1;
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
            if (!strcmp(part_type, "classic")) {
                METIS_PartMeshDual(&nelems, &nnodes, elem_ptr, elem_idx, NULL, NULL, &ncommon, &nparts, NULL, NULL, 
                        &objval, elem_part, node_part);
            } else {
                METIS_PartMeshNodal(&nelems, &nnodes, elem_ptr, elem_idx, NULL, NULL, &nparts, NULL, NULL, &objval,
                        elem_part, node_part);
            }
            for (i=0; i<nelems; i++) {
                (*metis_idx)[i] = (int) elem_part[i];
            }
            for (i=0; i<nprocs; i++) {
                intcell_per_proc[i] = 0;
            }
            //if our domain can't be divided in equal parts (breakets are very important!!!)
            extcell_per_proc[nprocs-1] = (nextcf_g-nextci_g+1)-(nprocs-1)
                    *((nextcf_g-nextci_g+1+(nprocs-1))/nprocs);
            //FIXME:move outside of this function
            fill_local_global_index(nprocs, *local_global_index_g, nelems, *metis_idx, intcell_per_proc);
        }
        //FIXME:move outside of this function
        count_ext_cells(nprocs, *local_global_index_g, nintci_g, nintcf_g, nextci_g, nextcf_g,
                lcc_g, *metis_idx, intcell_per_proc, extcell_per_proc);
    }
    return 0;
}


int allocate_local_variables(char* read_type, int myrank, int nprocs, int *nintci, int *nintcf, 
        int *nextci,int *nextcf, int ***lcc, double **bs, double **be, double **bn, double **bw, 
        double **bl, double **bh, double **bp, double **su, int* points_count, int*** points, 
        int** elems, int **local_global_index, int *intcell_per_proc, int *extcell_per_proc, 
        int *local_global_index_g, int points_count_g) {
    int i=0;
    MPI_Status status;
    //TODO: replace strcmp with definitions
    if (!strcmp(read_type, "oneread")) {
        // Before we allocate, we need to know how much memory to allocate same for all processes
        *nintci=0;
        if (myrank == 0) {
            *nintcf = intcell_per_proc[0]-1;
            *nextci = intcell_per_proc[0];
            *nextcf = *nintcf + extcell_per_proc[0];
            *points_count = points_count_g;
            for (i=1; i<nprocs; i++) {
                MPI_Send(&intcell_per_proc[i], 1, MPI_INT, i, NINTCF_SEND_INDEX, MPI_COMM_WORLD);
                MPI_Send(&extcell_per_proc[i], 1, MPI_INT, i, NEXTCF_SEND_INDEX, MPI_COMM_WORLD);
                MPI_Send(&points_count_g, 1, MPI_INT, i, POINTSCOUNT_SEND_INDEX, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(nintcf,1, MPI_INT, 0, NINTCF_SEND_INDEX, MPI_COMM_WORLD, &status);
            MPI_Recv(nextcf,1, MPI_INT, 0, NEXTCF_SEND_INDEX, MPI_COMM_WORLD, &status);
            MPI_Recv(points_count,1, MPI_INT, 0, POINTSCOUNT_SEND_INDEX, MPI_COMM_WORLD, &status);
            *nextci = *nintcf;
            // We need to subtract, because we got the number of cells(not the last index!)
            *nintcf = *nintcf-1;
            // We need to subtract, because we got the number of cells(not the last index!)
            *nextcf = *nintcf + *nextcf;
        }
    } else {
        *nintci = 0;
        *nintcf = intcell_per_proc[myrank]-1;
        *nextci = intcell_per_proc[myrank];
        *nextcf = *nintcf+extcell_per_proc[myrank];
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
    if ((*local_global_index = (int *) malloc(((*nintcf)+1)*sizeof(int))) == NULL) {
        fprintf(stderr, "malloc(local_global_index) failed\n");
        return -1;
    }
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


int send_or_read_data(char* read_type, int myrank, int nprocs, int nintci, int nintcf, int nextci, 
        int nextcf, int **lcc, double *bs, double *be, double *bn, double *bw, double *bl, 
        double *bh, double *bp, double *su, int points_count, int** points, int* elems, 
        int *local_global_index, int *intcell_per_proc, int *extcell_per_proc, int nintci_g, 
        int nintcf_g, int nextci_g, int nextcf_g, int **lcc_g, double **bs_g, double **be_g, 
        double **bn_g, double **bw_g, double **bl_g, double **bh_g, double **bp_g, double **su_g, 
        int points_count_g, int** points_g, int **elems_g, int *local_global_index_g) {
    int k=0, i=0;
    MPI_Status status;
    //TODO: replace strcmp with definitions
    if (!strcmp(read_type, "oneread")) {
        if (myrank == 0) {
            //TODO:check allocation status
            sort_data_by_local_global_index(nintci_g, nintcf_g, nextci_g, nextcf_g, &*bs_g, &*be_g, 
                    &*bn_g, &*bw_g, &*bl_g, &*bh_g, &*bp_g, &*su_g, &*elems_g, local_global_index_g);

            // Used to send data, from some index(because nodes can have not equal number of cells)
            int start_idx = 0;
            // Copy memory for process 0
            for(i=nintci; i<nintcf+1; i++) {
                memcpy(lcc[i], lcc_g[local_global_index_g[i]], 6*sizeof(int));
            }
            memcpy(bs, *bs_g, intcell_per_proc[0]*sizeof(double));
            memcpy(be, *be_g, intcell_per_proc[0]*sizeof(double));
            memcpy(bn, *bn_g, intcell_per_proc[0]*sizeof(double));
            memcpy(bw, *bw_g, intcell_per_proc[0]*sizeof(double));
            memcpy(bl, *bl_g, intcell_per_proc[0]*sizeof(double));
            memcpy(bh, *bh_g, intcell_per_proc[0]*sizeof(double));
            memcpy(bp, *bp_g, intcell_per_proc[0]*sizeof(double));
            memcpy(su, *su_g, intcell_per_proc[0]*sizeof(double));
            for (i=0; i<points_count; i++) {
                memcpy(points[i], points_g[i], 3*sizeof(int));
            }
            memcpy(elems, *elems_g, intcell_per_proc[0]*8*sizeof(int));
            memcpy(local_global_index, local_global_index_g, intcell_per_proc[0]*sizeof(int));
            // Send all other data
            for (k=1; k<nprocs; ++k) {
                start_idx += intcell_per_proc[k-1];
                for (i=0; i < intcell_per_proc[k]; ++i) {
                    MPI_Send(lcc_g[local_global_index_g[start_idx+i]], 6, MPI_INT, k, 0, MPI_COMM_WORLD);
                }
                MPI_Send(&(*bs_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
                MPI_Send(&(*be_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 2, MPI_COMM_WORLD);
                MPI_Send(&(*bn_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 3, MPI_COMM_WORLD);
                MPI_Send(&(*bw_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 4, MPI_COMM_WORLD);
                MPI_Send(&(*bl_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 5, MPI_COMM_WORLD);
                MPI_Send(&(*bh_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 6, MPI_COMM_WORLD);
                MPI_Send(&(*bp_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 7, MPI_COMM_WORLD);
                MPI_Send(&(*su_g)[start_idx], intcell_per_proc[k], MPI_DOUBLE, k, 8, MPI_COMM_WORLD);
                for (i=0; i<points_count; i++ ) {
                    MPI_Send(points_g[i], 3, MPI_INT, k, 99, MPI_COMM_WORLD);
                }
                MPI_Send(&(*elems_g)[start_idx*8], intcell_per_proc[k]*8, MPI_INT, k, 100, MPI_COMM_WORLD);
                MPI_Send(&local_global_index_g[start_idx], intcell_per_proc[k], MPI_INT, k, 101, MPI_COMM_WORLD);
            }
        } else {
            // Receive data from process 0
            for (i=nintci; i < nintcf; ++i) {
                MPI_Recv(lcc[i], 6, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            }
            MPI_Recv(bs, (nintcf-nintci+1), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(be, (nintcf-nintci+1), MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(bn, (nintcf-nintci+1), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
            MPI_Recv(bw, (nintcf-nintci+1), MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);
            MPI_Recv(bl, (nintcf-nintci+1), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status);
            MPI_Recv(bh, (nintcf-nintci+1), MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &status);
            MPI_Recv(bp, (nintcf-nintci+1), MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, &status);
            MPI_Recv(su, (nintcf-nintci+1), MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, &status);
            for (i=0; i<points_count; i++) {
                MPI_Recv(points[i], 3, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
            }
            MPI_Recv(elems, (nintcf-nintci+1)*8, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
            MPI_Recv(local_global_index, (nintcf-nintci+1), MPI_INT, 0, 101, MPI_COMM_WORLD, &status);
        }
    } else {
        //TODO: check allocation status
        sort_data_by_local_global_index(nintci_g, nintcf_g, nextci_g, nextcf_g, &*bs_g, &*be_g, 
                &*bn_g, &*bw_g, &*bl_g, &*bh_g, &*bp_g, &*su_g, &*elems_g, local_global_index_g);

        // Used to send data, from some index(because nodes can have not equal number of cells)
        int start_idx = 0;
        for (k=0; k<myrank; ++k) {
            start_idx += intcell_per_proc[k];
        }
        // Copy memory for process 0
        for(i=nintci; i<nintcf+1; i++) {
            memcpy(lcc[i], lcc_g[local_global_index_g[start_idx+i]], 6*sizeof(int));
        }
        memcpy(bs, &(*bs_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(be, &(*be_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(bn, &(*bn_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(bw, &(*bw_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(bl, &(*bl_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(bh, &(*bh_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(bp, &(*bp_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        memcpy(su, &(*su_g)[start_idx], intcell_per_proc[myrank]*sizeof(double));
        for (i=0; i<points_count; i++) {
            memcpy(points[i], points_g[i], 3*sizeof(int));
        }
        memcpy(elems, &(*elems_g)[start_idx*8], intcell_per_proc[myrank]*8*sizeof(int));
        memcpy(local_global_index, &local_global_index_g[start_idx], intcell_per_proc[myrank]*sizeof(int));
    }
    return 0;
}


/**
 * Compute local_global_index in sorted order
 */
//FIXME: change name of this function
void fill_local_global_index(int nprocs,int *local_global_index_g, int ne,  int *metis_idx, 
        int *intcell_per_proc) {
    int i=0, k=0, writer_counter=0;
    //TODO:split this for loop into 2 simpler ones
    for(k=0; k<nprocs; ++k) {
        for(i=0; i<ne; ++i) {
            if(metis_idx[i]==k) {
                ++intcell_per_proc[k];
                local_global_index_g[writer_counter]=i;
                ++writer_counter;
            }
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


void count_ext_cells(int nprocs, int *local_global_index_g, int nintci_g, int nintcf_g, int nextci_g, 
        int nextcf_g, int **lcc_g, int *metis_idx, int *intcell_per_proc, int *extcell_per_proc) {
    int is_int_cell=0, n_ghost_cells=0, start_idx=0, idx=0, i=0, j=0, proc=0;
    int g2l[nprocs][nextcf_g+1], isSaved[nextcf_g+1];
    
    for (i=0; i<nextcf_g+1; ++i) {
        for (proc=0; proc<nprocs; ++proc) {
            g2l[proc][i]=-1;
        }
    }

    for (proc=0; proc<nprocs; ++proc) {
    	if (proc != 0) {
            start_idx += intcell_per_proc[proc-1];
        }
    	// Set all isUsed to zero
    	memset(isSaved, 0, (nextcf_g+1)*sizeof(int));
    	// Compute external cell
    	extcell_per_proc[proc] = 0;	// TODO: is it good place for it?
    	for (i=0; i<intcell_per_proc[proc]; ++i) {
            idx = local_global_index_g[start_idx+i];
            for (j=0; j<6; ++j) {
                if (lcc_g[idx][j] >= nintcf_g+1 && !isSaved[lcc_g[idx][j]]) {
                    isSaved[lcc_g[idx][i]] = 1;
                    lcc_g[idx][j] = g2l[proc][idx] = intcell_per_proc[proc]+extcell_per_proc[proc];
                    ++extcell_per_proc[proc];
                }
            }
    	}
    	// End compute external cell
    	memset(isSaved, 0, (nextcf_g+1)*sizeof(int));
    	// Ghost cell and local indexing
    	n_ghost_cells = 0;
    	for (i=0; i<intcell_per_proc[proc]; ++i) {
            idx = local_global_index_g[start_idx + i];
            for (j=0; j<6; ++j) {
                // Check if this cell is in this processor
                is_int_cell = lcc_g[idx][j] <= nintcf_g && metis_idx[lcc_g[idx][j]] == proc;
                if (!is_int_cell && !isSaved[lcc_g[idx][j]]) {
                    isSaved[lcc_g[idx][i]] = 1;
                    lcc_g[idx][j] = g2l[proc][idx] = intcell_per_proc[proc]+extcell_per_proc[proc]
                            +n_ghost_cells;
                    ++n_ghost_cells;
                } else {
                    lcc_g[idx][j] = g2l[proc][idx] = i;
                }
            }
    	}
    	// End ghost cell
    	extcell_per_proc[proc] += n_ghost_cells;
    }
}
