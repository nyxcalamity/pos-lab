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
int allocate_lcc_elems_points(int read_key, int myrank, int nprocs, int *nintci, int *nintcf, 
        int ***lcc, int* points_count, int*** points, int** elems, int **local_global_index, 
        int points_count_g);


/**
 * Exchange or copy from global all local data.
 */
int fill_lcc_elems_points(int read_key, int myrank, int nprocs, int nintci, int nintcf, int **lcc, 
        int points_count, int** points, int* elems, int *local_global_index,  int **lcc_g, 
        int points_count_g, int** points_g, int **elems_g);

/**
 * Allocates local memory for various boundary coefficients.
 */
int allocate_boundary_coef(int read_key, int myrank, int nprocs, int *nextcf, double **bs, 
        double **be, double **bn, double **bw, double **bl, double **bh, double **bp, double **su);


/**
 * Exchange local values of boundary coefficients or copy them from global values.
 */
int fill_boundary_coef(int read_key, int myrank, int nprocs, int nintci, int nintcf, int nextci,
        int nextcf, double *bs, double *be, double *bn, double *bw, double *bl,  double *bh, 
        double *bp, double *su, int *local_global_index, double **bs_g, double **be_g, double **bn_g, 
        double **bw_g, double **bl_g, double **bh_g, double **bp_g, double **su_g);


/**
 * Computes a 1D array with its indexes as local cell IDs and its values as global cell IDs.
 */
void fill_l2g(int myrank, int nintcf,  int** local_global_index, int *partitioning_map, int nelems_g);


/**
 * Counts the number of external cell for each process and change lcc to have indexes which are local
 *
 * FIXME: delete unused arguments + improve comments
 * FIXME: check on bad allocation
 */
void build_lists_g2l_next(int nprocs, int myrank, int *partitioning_map, 
        int nintcf_g, int nextci_g, int nextcf_g, int* nintci, int* nintcf, int* nextci, 
        int* nextcf, int*** lcc, int* points_count, int*** points, int** elems, double** var, 
        double** cgup, double** oc, double** cnorm, int** local_global_index, int** global_local_index, 
        int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, int **recv_cnt, int*** recv_lst);


/**
 * Sends and receives sizes of send list(which is equal to size of receive list) and allocate needed memory
 *
 * FIXME: improve comments
 */
void allocate_send_lists(int myrank, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
        int **recv_cnt);


/**
 * Send receive lists and get them in send lists.
 *
 * FIXME: better comment + delete unused arguments
 */
void exchange_lists(int myrank, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
        int **recv_cnt, int*** recv_lst);

#endif /* INITIALIZATION_ALGORITHMS_H_ */
