#ifndef INITIALIZATION_ALGORITHMS_H_
#define INITIALIZATION_ALGORITHMS_H_


/**
 * Reads data needed for computation.
 */
int read_init_data(char* file_in, int read_key, int myrank, int *nintci, int *nintcf, 
        int *nextci, int *nextcf, int ***lcc, double **bs, double **be, double **bn, double **bw, 
        double **bl, double **bh, double **bp, double **su, int* points_count, int***points, int** elems);


/**
 * Computes external and internal cells starting and ending positions.
 */
int partition(int part_key, int read_key, int myrank, int nprocs, int nintci_g,
        int nintcf_g, int nextci_g, int nextcf_g, int *nintci, int *nintcf, int *nextci, int *nextcf,
        int **lcc_g, int points_count_g, int**points_g, int* elems_g, int *intcell_per_proc,
        int *extcell_per_proc, int** local_global_index_g, int** local_global_index, int **partitioning_map);


/**
 * Exchange or find in global array nintci,nintcf,nextci,nextcf and allocate corresponding memory for
 * elements, points and lcc.
 */
int allocate_lcc_elems_points(int read_key, int myrank, int nprocs, int *nintci, int *nintcf, int *nextci,
        int ***lcc, int* points_count, int*** points, int** elems, int **local_global_index, 
        int points_count_g, int *int_cells_per_proc);


/**
 * Exchange or copy from global all local data.
 */
int fill_lcc_elems_points(int read_key, int myrank, int nprocs, int nintci, int nintcf, int **lcc, 
        int points_count, int** points, int* elems, int *local_global_index, int **local_global_index_g, 
        int **lcc_g, int points_count_g, int** points_g, int **elems_g, int *int_cells_per_proc);

/**
 * Allocates local memory for various boundary coefficients.
 */
int allocate_boundary_coef(int *nextcf, double **bs, double **be, double **bn, double **bw, double **bl, 
        double **bh, double **bp, double **su);


/**
 * Exchange local values of boundary coefficients or copy them from global values.
 */
int fill_boundary_coef(int read_key, int myrank, int nprocs, int nintci, int nintcf, int nextci,
        int nextcf, double *bs, double *be, double *bn, double *bw, double *bl,  double *bh, 
        double *bp, double *su, int *local_global_index, int **local_global_index_g, double **bs_g, 
        double **be_g, double **bn_g, double **bw_g, double **bl_g, double **bh_g, double **bp_g, 
        double **su_g, int *int_cells_per_proc);


/**
 * Computes a 1D array with its indexes as local cell IDs and its values as global cell IDs.
 */
void fill_l2g(int read_key, int myrank, int nproc, int nintcf,  int** local_global_index, 
        int ***local_global_index_g, int *partitioning_map, int nelems_g, int *int_cells_per_proc);


/**
 * Counts the number of external cell for each process and change lcc to have indexes which are local
 */
void build_lists_g2l_next(int nprocs, int myrank, int *partitioning_map, int nintcf_g, int nextcf_g, 
        int* nintcf, int* nextcf, int*** lcc, int** local_global_index, int** global_local_index, 
        int *nghb_cnt, int** nghb_to_rank, int **recv_cnt, int*** recv_lst);


/**
 * Sends and receives sizes of send list(which is equal to size of receive list) and allocate needed memory
 */
void allocate_send_lists(int myrank, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
        int **recv_cnt);


/**
 * Send receive lists and get them in send lists.
 */
void exchange_lists(int myrank, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
        int **recv_cnt, int*** recv_lst);

/**
 * Broadcasts partitioning map to all processors, if necessary.
 */
void bcast_partitioning(int read_key, int myrank, int **partitioning_map, int *nintci_g, int *nintcf_g,
        int *nextci_g, int *nextcf_g);


/**
 * Parses lcc and communication lists in order to replace global cells indexes to local ones.
 */
void converte_global2local_idx(int myrank, int *g2l, int nintci, int nintcf, int **lcc, 
        int ngbh_cnt, int *send_cnt, int **send_lst, int *recv_cnt, int **recv_lst);

#endif /* INITIALIZATION_ALGORITHMS_H_ */
