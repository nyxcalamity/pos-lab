#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "metis.h"

#include "initialization_algorithms.h"
#include "util_read_files.h"
#include "posl_definitions.h"


int read_init_data(char* file_in, int read_key, int myrank, int *nintci, int *nintcf, 
        int *nextci, int *nextcf, int ***lcc, double **bs, double **be, double **bn, double **bw, 
        double **bl, double **bh, double **bp, double **su, int* points_count, int***points, 
        int** elems) {
    int f_status=0;
    if (read_key == POSL_INIT_ONE_READ) {
        if (myrank == 0) {
            f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
                    &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count, &*points, &*elems);
        }
    } else {
        //TODO:implement sheer geometry reading
        f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs, 
                &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count, &*points, &*elems);
    }
    
    return f_status;
}


void bcast_partitioning(int read_key, int myrank, int **partitioning_map, int *nintci_g, int *nintcf_g,
        int *nextci_g, int *nextcf_g){
    if (read_key == POSL_INIT_ONE_READ) {
        //broadcast partitioning map size
        //TODO:this is same as broadcasting cell idx sizes, merge the two
        if (myrank != 0) {
            *nintci_g = 0; *nintcf_g = 0;
            *nextci_g = 0; *nextcf_g = 0;
        }
        MPI_Bcast(nintci_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(nintcf_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(nextci_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(nextcf_g, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        int partitioning_size = (*nintcf_g-*nintci_g+1);
        
        //broadcast partition map itself
        if (myrank != 0) {
            *partitioning_map = (int *) calloc(partitioning_size, sizeof(int));
        }
        MPI_Bcast(*partitioning_map, partitioning_size, MPI_INT, 0, MPI_COMM_WORLD);
    }
}


int partition(int part_key, int read_key, int myrank, int nprocs, int nintci_g, 
        int nintcf_g, int nextci_g, int nextcf_g, int *nintci, int *nintcf, int *nextci, int *nextcf, 
        int **lcc_g, int points_count_g, int **points_g, int *elems_g, int *int_cells_per_proc, 
        int *extcell_per_proc, int **local_global_index_g, int **local_global_index, int **partitioning_map) {
    int i=0;
    idx_t nelems, nnodes, ncommon, nparts, objval;
    idx_t *elem_ptr, *elem_idx, *elem_part, *node_part;
    
    nelems = nintcf_g-nintci_g+1;
    *partitioning_map = (int *) calloc(sizeof(int), (nintcf_g-nintci_g+1));
    
    if (((read_key == POSL_INIT_ONE_READ) && (myrank == 0)) || (read_key == POSL_INIT_ALL_READ)) {
        *nintci = 0; *nintcf = 0;
        if (part_key == POSL_PARTITIONING_CLASSIC) {
            //the last processor always gets different number of cells
            int elem_per_proc = (nelems+(nprocs-1))/nprocs;
            *nextci = (myrank == nprocs-1) ? nelems-(nprocs-1)*elem_per_proc : elem_per_proc;
            *nintcf = *nextci-1;
            
            //build global cell allocation
            if (read_key == POSL_INIT_ONE_READ) {
                for (i=0; i<(nprocs-1); ++i) {
                    int_cells_per_proc[i] = elem_per_proc;
                }
                int_cells_per_proc[nprocs-1] = nelems-(nprocs-1)*elem_per_proc;
            }
                        
            for (i=0; i<nelems; ++i) {
                (*partitioning_map)[i] = i/elem_per_proc;
            }
        } else {
            //initialize variables for metis
            nnodes = points_count_g;
            //FIXME:perform a proper constant externalization
            ncommon = 4;
            nparts = nprocs;
            elem_ptr = (idx_t *) calloc(nelems+1, sizeof(idx_t));
            elem_idx = (idx_t *) calloc(nelems*8, sizeof(idx_t));
            elem_part = (idx_t *) calloc(nelems, sizeof(idx_t));
            node_part = (idx_t *) calloc(nnodes, sizeof(idx_t));
            //TODO:add comments regarding element allocation and input data for metis
            for (i=0; i<(nelems+1); i++) {
                elem_ptr[i] = 8*i;
            }
            for (i=0; i<(nelems*8); i++) {
                elem_idx[i] = elems_g[i];
            }

            //perform metis partitioning
            if (part_key == POSL_PARTITIONING_DUAL) {
                METIS_PartMeshDual(&nelems, &nnodes, elem_ptr, elem_idx, NULL, NULL, &ncommon, 
                        &nparts, NULL, NULL, &objval, elem_part, node_part);
            } else {
                METIS_PartMeshNodal(&nelems, &nnodes, elem_ptr, elem_idx, NULL, NULL, &nparts, 
                        NULL, NULL, &objval, elem_part, node_part);
            }
            
            for (i=0; i<nelems; i++) {
                (*partitioning_map)[i] = (int) elem_part[i];
            }
            
            //initialize global cell counters
            if (read_key == POSL_INIT_ONE_READ) {
                for (i=0; i<nprocs; ++i) {
                    int_cells_per_proc[i] = 0;
                }
            }
            
            //compute position of last internal cell
            for (i=0; i<nelems; i++) {
                if (read_key == POSL_INIT_ONE_READ) {
                    int_cells_per_proc[(*partitioning_map)[i]] += 1;
                } else {
                    if (myrank == (*partitioning_map)[i]) {
                        (*nintcf) += 1;
                    }
                }
            }
            
            //assign local internal cell ending idx
            if (read_key == POSL_INIT_ONE_READ) {
                *nintcf = int_cells_per_proc[myrank];
            }
            *nextci = (*nintcf)--;
        }
    }
    return 0;
}


int allocate_lcc_elems_points(int read_key, int myrank, int nprocs, int *nintci, int *nintcf, 
        int ***lcc, int* points_count, int*** points, int** elems, int **local_global_index, 
        int points_count_g, int *int_cells_per_proc) {
    int i=0;
    MPI_Status status;
    //starting index is same for all processors
    *nintci = 0;
    if (read_key == POSL_INIT_ONE_READ) {
        if (myrank == 0) {
            *points_count = points_count_g;
            for (i=1; i<nprocs; i++) {
                MPI_Send(&int_cells_per_proc[i], 1, MPI_INT, i, POSL_MPI_TAG_NINTCF, MPI_COMM_WORLD);
                MPI_Send(&points_count_g, 1, MPI_INT, i, POSL_MPI_TAG_POINTS_COUNT, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(nintcf, 1, MPI_INT, 0, POSL_MPI_TAG_NINTCF, MPI_COMM_WORLD, &status);
            MPI_Recv(points_count, 1, MPI_INT, 0, POSL_MPI_TAG_POINTS_COUNT, MPI_COMM_WORLD, &status);
            //decrement the value, since our indexing starts from zero
            --(*nintcf);
        }
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
    if ((*points = (int **) calloc(*points_count, sizeof(int*))) == NULL) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }
    for (i=0; i<*points_count; i++) {
        if (((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }
    
    return 0;
}


int fill_lcc_elems_points(int read_key, int myrank, int nprocs, int nintci, int nintcf, int **lcc, 
        int points_count, int** points, int* elems, int *local_global_index,  int **local_global_index_g,
        int **lcc_g, int points_count_g, int** points_g, int **elems_g, int *int_cells_per_proc) {
    int proc=0, i=0;
    int *buff;
    MPI_Status status;
    if (read_key == POSL_INIT_ONE_READ) {
        if (myrank == 0) {
            for (proc=1; proc<nprocs; ++proc) {
                for (i=0; i<int_cells_per_proc[proc]; ++i) {
                    MPI_Send(lcc_g[local_global_index_g[proc][i]], 6, MPI_INT, proc, 
                            POSL_MPI_TAG_LCC, MPI_COMM_WORLD);
                }
                
                for (i=0; i<points_count; i++) {
                    MPI_Send(points_g[i], 3, MPI_INT, proc, POSL_MPI_TAG_POINTS, MPI_COMM_WORLD);
                }
                
                //TODO:Use MPI custom data types (indexing))
                //allocate proper buffer memory
                if ((buff = (int *) malloc(int_cells_per_proc[proc]*8*sizeof(int))) == NULL) {
                    fprintf(stderr, "malloc() element buffer allocation failed\n");
                    return -1;
                }
                
                //fill the buffer
                for (i=0; i<int_cells_per_proc[proc]; ++i) {
                    memcpy(&(buff[8*i]), &((*elems_g)[8*local_global_index_g[proc][i]]), 8*sizeof(int));
                }
                                
                //send the elements
                MPI_Send(buff, int_cells_per_proc[proc]*8, MPI_INT, proc, POSL_MPI_TAG_ELEMENTS, 
                        MPI_COMM_WORLD);
                
                //free buffer memory
                free(buff);
            }
        } else {
            for (i=nintci; i<nintcf; ++i) {
                MPI_Recv(lcc[i], 6, MPI_INT, 0, POSL_MPI_TAG_LCC, MPI_COMM_WORLD, &status);
            }
            for (i=0; i<points_count; i++) {
                MPI_Recv(points[i], 3, MPI_INT, 0, POSL_MPI_TAG_POINTS, MPI_COMM_WORLD, &status);
            }
            MPI_Recv(elems, (nintcf+1)*8, MPI_INT, 0, POSL_MPI_TAG_ELEMENTS, MPI_COMM_WORLD, &status);
        }
    }
    if (read_key == POSL_INIT_ALL_READ || (read_key == POSL_INIT_ONE_READ && myrank == 0)) {
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


int allocate_boundary_coef(int *nextcf, double **bs, double **be, double **bn, double **bw, double **bl, 
        double **bh, double **bp, double **su) {
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


int fill_boundary_coef(int read_key, int myrank, int nprocs, int nintci, int nintcf, int nextci,
        int nextcf, double *bs, double *be, double *bn, double *bw, double *bl,  double *bh, 
        double *bp, double *su, int *local_global_index, int **local_global_index_g, double **bs_g, 
        double **be_g, double **bn_g, double **bw_g, double **bl_g, double **bh_g, double **bp_g, 
        double **su_g, int *int_cells_per_proc) {
    int proc=0, i=0;
    double *buf_bs, *buf_be, *buf_bn, *buf_bw, *buf_bl, *buf_bh, *buf_bp, *buf_su;
    MPI_Status status;
    if (read_key == POSL_INIT_ONE_READ) {
        if (myrank == 0) {
            for (proc=1; proc<nprocs; ++proc) {
                //allocate buffers
                if ((buf_bs = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_be = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_bn = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_bw = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_bl = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_bh = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_bp = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }
                if ((buf_su = (double *) malloc(int_cells_per_proc[proc]*sizeof(double))) == NULL) {
                    printf("malloc() b value failed\n");
                    return -1;
                }

                //fill in buffer data
                for (i=0; i<int_cells_per_proc[proc]; ++i) {
                    memcpy(&(buf_bs[i]), &(*bs_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_be[i]), &(*be_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_bn[i]), &(*bn_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_bw[i]), &(*bw_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_bl[i]), &(*bl_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_bh[i]), &(*bh_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_bp[i]), &(*bp_g)[local_global_index_g[proc][i]], sizeof(double));
                    memcpy(&(buf_su[i]), &(*su_g)[local_global_index_g[proc][i]], sizeof(double));
                }
                
                //send the data
                //FIXME: externalize mpi tags
                MPI_Send(buf_bs, int_cells_per_proc[proc], MPI_DOUBLE, proc, 101, MPI_COMM_WORLD);
                MPI_Send(buf_be, int_cells_per_proc[proc], MPI_DOUBLE, proc, 102, MPI_COMM_WORLD);
                MPI_Send(buf_bn, int_cells_per_proc[proc], MPI_DOUBLE, proc, 103, MPI_COMM_WORLD);
                MPI_Send(buf_bw, int_cells_per_proc[proc], MPI_DOUBLE, proc, 104, MPI_COMM_WORLD);
                MPI_Send(buf_bl, int_cells_per_proc[proc], MPI_DOUBLE, proc, 105, MPI_COMM_WORLD);
                MPI_Send(buf_bh, int_cells_per_proc[proc], MPI_DOUBLE, proc, 106, MPI_COMM_WORLD);
                MPI_Send(buf_bp, int_cells_per_proc[proc], MPI_DOUBLE, proc, 107, MPI_COMM_WORLD);
                MPI_Send(buf_su, int_cells_per_proc[proc], MPI_DOUBLE, proc, 108, MPI_COMM_WORLD);
                
                //free buffer memory
                free(buf_bs);
                free(buf_bn);
                free(buf_bw);
                free(buf_bl);
                free(buf_bh);
                free(buf_bp);
                free(buf_su);
            }
        } else {
            MPI_Recv(bs, (nintcf+1), MPI_DOUBLE, 0, 101, MPI_COMM_WORLD, &status);
            MPI_Recv(be, (nintcf+1), MPI_DOUBLE, 0, 102, MPI_COMM_WORLD, &status);
            MPI_Recv(bn, (nintcf+1), MPI_DOUBLE, 0, 103, MPI_COMM_WORLD, &status);
            MPI_Recv(bw, (nintcf+1), MPI_DOUBLE, 0, 104, MPI_COMM_WORLD, &status);
            MPI_Recv(bl, (nintcf+1), MPI_DOUBLE, 0, 105, MPI_COMM_WORLD, &status);
            MPI_Recv(bh, (nintcf+1), MPI_DOUBLE, 0, 106, MPI_COMM_WORLD, &status);
            MPI_Recv(bp, (nintcf+1), MPI_DOUBLE, 0, 107, MPI_COMM_WORLD, &status);
            MPI_Recv(su, (nintcf+1), MPI_DOUBLE, 0, 108, MPI_COMM_WORLD, &status);
        }
    }
    if (read_key == POSL_INIT_ALL_READ || (read_key == POSL_INIT_ONE_READ && myrank ==0)) {
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


void fill_l2g(int read_key, int myrank, int nproc, int nintcf, int** local_global_index, 
        int ***local_global_index_g, int *partitioning_map, int nelems_g, int *int_cells_per_proc) {
    int i=0, local_idx=0, current_proc=0;
    if ((*local_global_index = (int *) malloc(((nintcf)+1)*sizeof(int))) == NULL) {
        //FIXME:handle this error
        fprintf(stderr, "malloc(local_global_index) failed\n");
    }

    for (i=0; i<nelems_g; ++i) {
        if (partitioning_map[i] == myrank) {
            (*local_global_index)[local_idx] = i;
            ++local_idx;
        }
    }
    
    //TODO:refactor code to reduce duplicates
    //in addition build global global to local if applicable
    if (read_key == POSL_INIT_ONE_READ && myrank == 0) {        
        //FIXME:handle memory allocation errors
        if ((*local_global_index_g = (int**) malloc(nproc*sizeof(int*))) == NULL){
            fprintf(stderr, "malloc failed to allocate first dimension of l2g_g");
        }
        for (i=0; i<nproc; i++) {
            if (((*local_global_index_g)[i] = (int*) malloc(int_cells_per_proc[i]*sizeof(int))) == NULL){
                fprintf(stderr, "malloc failed to allocate second dimension of l2g_g");
            }
        }

        //initialize array of local indexes
        int last_local_idx[nproc];
        for (i=0; i<nproc; ++i) {
            last_local_idx[i] = 0;
        }
        
        for (i=0; i<nelems_g; ++i) {
            current_proc = partitioning_map[i];
            local_idx = last_local_idx[current_proc];
            
            if (local_idx > int_cells_per_proc[current_proc]) {
                printf("Wrong local index for processor #%d", current_proc);
            }
            
            (*local_global_index_g)[current_proc][local_idx] = i;
            ++(last_local_idx[current_proc]);
        }
    }
}


void build_lists_g2l_next(int nprocs, int myrank, int *partitioning_map, int nintcf_g, int nextcf_g, 
        int* nintcf, int* nextcf, int*** lcc, int** local_global_index, int** global_local_index, 
        int *nghb_cnt, int** nghb_to_rank, int **recv_cnt, int*** recv_lst) {
    /*********** Initialize and allocate g2l and tmp variables ************************************/
    int proc=0, i=0, j=0, is_ghost_cell=0, idx_g=-1, neighbor_rank=-1; /// -1 to track errors
    int n_ghost_cells[nprocs], start_idx_per_proc[nprocs];
    int nghb_idx=0;    /// Needed for filling of nghb_to_rank
    memset(n_ghost_cells, 0, nprocs*sizeof(int));
    memset(start_idx_per_proc, 0, nprocs*sizeof(int));
    *nextcf = 0;
    // Used to check if we already saved the cel index
    int is_saved[nextcf_g+1];
    /** Lists for neighboring information */
    /// number of cells to be received from each neighbour (size: nprocs)
    int tmp_recv_cnt[nprocs];
    memset(tmp_recv_cnt, 0, nprocs*sizeof(int));
    // TODO: try not to use tmp_recv_lst and save all data immediately in revc_lst
    /// lists of cells to be received from each neighbour (size: nprocs x recv_cnt[*])
    int *tmp_recv_lst[nprocs];

    if ((*global_local_index = (int *) malloc((nextcf_g+1)*sizeof(int))) == NULL) {
        fprintf(stderr, "malloc(global_local_index) in build_lists_g2l_next failed\n");
        //FIXME:handle this error
    }
    // We will fill whole g2l by -1 for easier debugging
    for (i=0; i<=nextcf_g; ++i) {
        (*global_local_index)[i] = -1;
    }
    /*********** End initialize and allocate g2l and tmp variables ********************************/
// TODO: We can do all work in one big loop but then we need to use the buffer before we allocate
//    memory for tmp_recv_lst
    /************************ Start processing lcc and generating data ****************************/
    memset(is_saved, 0, (nextcf_g+1)*sizeof(int));
    // Count number of external cells and fill g2l with external cells, and fill tmp_recv_cnt
    for (i=0; i<=(*nintcf); ++i) {
        for (j=0; j<6; ++j) {
            idx_g = (*lcc)[i][j];
            if(idx_g<=nintcf_g) {
                neighbor_rank = partitioning_map[idx_g];
                is_ghost_cell = neighbor_rank != myrank;
            }

            if (!is_saved[idx_g]) {
                if((*lcc)[i][j] > nintcf_g) {
                    ++(*nextcf);
                    is_saved[idx_g] = 1;
                    (*global_local_index)[idx_g] = (*nintcf) + (*nextcf);
                } else if (is_ghost_cell) {
                    ++tmp_recv_cnt[neighbor_rank];
                    is_saved[idx_g] = 1;
                }
            }
        }
        // Fill g2l with internal cell
        (*global_local_index)[(*local_global_index)[i]] = i;
    }
    /// After the loop nextcf shows amount of elements
    *nextcf = *nintcf + *nextcf;

    // Allocate tmp_recv_lst
    for (proc=0; proc<nprocs; ++proc) {
        if(tmp_recv_cnt[proc] == 0) {
            tmp_recv_lst[proc] = NULL;
        } else {
            tmp_recv_lst[proc] = (int *) calloc(sizeof(int), tmp_recv_cnt[proc]);
        }
    }
    // Fill tmp_recv_list and g2l with ghost cells
    start_idx_per_proc[0]=0;
    for(proc=1; proc<nprocs; ++proc) {
        start_idx_per_proc[proc] = start_idx_per_proc[proc-1] + tmp_recv_cnt[proc-1];
    }
    memset(is_saved, 0, (nextcf_g+1)*sizeof(int));
    for (i=0; i<=(*nintcf); ++i) {
        for (j=0; j<6; ++j) {
            idx_g = (*lcc)[i][j];
            if(idx_g<=nintcf_g) {
                neighbor_rank = partitioning_map[idx_g];
                is_ghost_cell = neighbor_rank != myrank;
            } else {
                is_ghost_cell = 0;
            }
            if (is_ghost_cell && !is_saved[idx_g]) {
                is_saved[idx_g] = 1;
                tmp_recv_lst[neighbor_rank][n_ghost_cells[neighbor_rank]] = idx_g;
                ++n_ghost_cells[neighbor_rank];
                (*global_local_index)[idx_g] = (*nextcf) +
                        start_idx_per_proc[neighbor_rank] + n_ghost_cells[neighbor_rank];
            }
        }
    }

    // Count neighbors it is already set to zero in gccg.c
    for (proc=0; proc<nprocs; ++proc) {
        if (tmp_recv_cnt[proc] !=0) ++(*nghb_cnt);
    }

    // Allocate and fill nhb_to_rank
    if ((*nghb_to_rank = (int *) malloc((*nghb_cnt)*sizeof(int))) == NULL) {
        fprintf(stderr, "malloc(nghb_to_rank) in build_lists_g2l_next failed\n");
        //FIXME:handle this error
    }
    for (proc=0;proc<nprocs; ++proc) {
        if(tmp_recv_cnt[proc] !=0) {
            (*nghb_to_rank)[nghb_idx] = proc;
            ++nghb_idx;
        }
    }
    // Allocate and fill recv_cnt and recv_lst with global lcc indexing for further use
    // FIXME: is this good allocation?
    *recv_cnt = (int*) calloc(sizeof(int), *nghb_cnt);
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        (*recv_cnt)[nghb_idx] = tmp_recv_cnt[ (*nghb_to_rank)[nghb_idx] ];
    }
    *recv_lst = (int**) malloc( (*nghb_cnt)*sizeof(int*) );
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        (*recv_lst)[nghb_idx] = (int *) malloc( (*recv_cnt)[nghb_idx] * sizeof(int) );
    }
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        memcpy( (*recv_lst)[nghb_idx], tmp_recv_lst[ (*nghb_to_rank)[nghb_idx] ],
                (*recv_cnt)[nghb_idx]*sizeof(int) );
    }

//    printf("r%d, recv_lst[0][0]=%d, recv_lst[0][1]=%d, recv_lst[0][2]=%d, recv_lst[0][3]=%d\n",
//            myrank, (*recv_lst)[0][0], (*recv_lst)[0][1], (*recv_lst)[0][2], (*recv_lst)[0][3]);
//    printf("r%d, recv_cnt[0]=%d, recv_cnt[1]=%d, recv_cnt[2]=%d, recv_cnt[3]=%d\n",
//            myrank, (*recv_cnt)[0], (*recv_cnt)[1], (*recv_cnt)[2], (*recv_cnt)[3]);
//    printf("r%d, tmp_recv_cnt[0]=%d, tmp_recv_cnt[1]=%d, tmp_recv_cnt[2]=%d, tmp_recv_cnt[3]=%d\n",
//            myrank, tmp_recv_cnt[0], tmp_recv_cnt[1], tmp_recv_cnt[2], tmp_recv_cnt[3]);
//    printf("r%d, n_ghost_cells[0]=%d, n_ghost_cells[1]=%d, n_ghost_cells[2]=%d, n_ghost_cells[3]=%d\n",
//            myrank, n_ghost_cells[0], n_ghost_cells[1], n_ghost_cells[2], n_ghost_cells[3]);
//    printf("r%d, start_idx_per_proc[0]=%d, start_idx_per_proc[1]=%d, start_idx_per_proc[2]=%d, start_idx_per_proc[3]=%d\n",
//            myrank, start_idx_per_proc[0], start_idx_per_proc[1], start_idx_per_proc[2], start_idx_per_proc[3]);
//TODO:    printf("r%d, nghb_cnt=%d\n", myrank, (*nghb_cnt));
//TODO:    printf("r%d, nghb_to_rank[0]=%d, nghb_to_rank[1]=%d, nghb_to_rank[2]=%d, nghb_to_rank[3]=%d, nghb_to_rank[4]=%d\n",
//            myrank, (*nghb_to_rank)[0], (*nghb_to_rank)[1], (*nghb_to_rank)[2], (*nghb_to_rank)[3], (*nghb_to_rank)[3]);

    // Total number of external cells
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        *nextcf += (*recv_cnt)[nghb_idx];
    }
    /************************ End processing lcc and generating data ******************************/
    // Free memory
    for (proc=0; proc<nprocs; ++proc) {
        free(tmp_recv_lst[proc]);
    }
}


void allocate_send_lists(int myrank, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, 
        int*** send_lst, int **recv_cnt) {
    MPI_Status status;
    int nghb_idx=0;
    // Allocate send_cnt
    *send_cnt = (int*) calloc(sizeof(int), *nghb_cnt);
    // Send size of receive list for that processor, it will be equal to the size of send list of
    // received processor
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        // Send sizes
        MPI_Send(&(*recv_cnt)[nghb_idx], 1, MPI_INT, (*nghb_to_rank)[nghb_idx], myrank, MPI_COMM_WORLD);
        // Receive sizes
        MPI_Recv(&(*send_cnt)[nghb_idx],1 , MPI_INT, (*nghb_to_rank)[nghb_idx], 
                (*nghb_to_rank)[nghb_idx], MPI_COMM_WORLD, &status);
    }
    // Allocate send_lst with given sizes
    *send_lst = (int**) malloc( (*nghb_cnt)*sizeof(int*) );
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        (*send_lst)[nghb_idx] = (int *) malloc((*send_cnt)[nghb_idx]*sizeof(int));
    }
}


void exchange_lists(int myrank, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
        int **recv_cnt, int*** recv_lst) {
    MPI_Status status;
    int nghb_idx=0;
    MPI_Request request_send[*nghb_cnt], request_recv[*nghb_cnt];
    // Send receive list and save it in send list
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        // Send sizes
        MPI_Isend((*recv_lst)[nghb_idx], (*recv_cnt)[nghb_idx], MPI_INT, (*nghb_to_rank)[nghb_idx],
                myrank, MPI_COMM_WORLD, &request_send[nghb_idx]);
        // Receive sizes
        MPI_Irecv((*send_lst)[nghb_idx], (*send_cnt)[nghb_idx], MPI_INT, (*nghb_to_rank)[nghb_idx], 
                (*nghb_to_rank)[nghb_idx], MPI_COMM_WORLD, &request_recv[nghb_idx]);
    }
    // Synchronize everything
    for (nghb_idx=0; nghb_idx<(*nghb_cnt); ++nghb_idx) {
        MPI_Wait(&request_send[nghb_idx], &status);
        MPI_Wait(&request_recv[nghb_idx], &status);
    }
}
