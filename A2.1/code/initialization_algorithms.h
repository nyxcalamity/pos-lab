#ifndef INITIALIZATION_ALGORITHMS_H_
#define INITIALIZATION_ALGORITHMS_H_

/**
 * Reads global data if we have "oneread" read type or
 * reads data needed by METIS computation if we have "allread"
 *
 * @param file_name
 * @param read_type
 * @param myrank
 * @param nintci
 * @param nintcf
 * @param nextci
 * @param nextcf
 * @param lcc
 * @param BS
 * @param BE
 * @param BN
 * @param BW
 * @param BL
 * @param BH
 * @param BP
 * @param SU
 * @param points_count
 * @param points
 * @param elems
 * @return
 */
int read_global_data_or_geometry(char* file_in,char* read_type, int myrank,
        int *nintci, int *nintcf, int *nextci, int *nextcf,
        int ***lcc,
        double **bs, double **be, double **bn, double **bw,
		double **bl, double **bh, double **bp, double **su,
		int* points_count, int***points, int** elems);

/**
 * Uses read geometry data to allocate and compute local_global_index and
 * to fill cell_per_node
 *
 * @param file_name
 * @param read_type
 * @param myrank
 * @param nprocs
 * @param nintci_g
 * @param nintcf_g
 * @param nextci_g
 * @param nextcf_g
 * @param nintci
 * @param nintcf
 * @param nextci
 * @param nextcf
 * @param lcc
 * @param points_count
 * @param points
 * @param elems
 * @param intcell_per_proc
 * @param extcell_per_proc
 * @param local_global_index
 * @return
 */
int compute_metis(char* part_type, char* read_type, int myrank, int nprocs,
        int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
		int *nintci, int *nintcf, int *nextci, int *nextcf,
        int **lcc_g,
		int points_count_g, int**points_g, int* elems_g,
		int *intcell_per_proc, int *extcell_per_proc, int** local_global_index,
		int **metis_idx);

/**
 * Send and receive or compute nintci,nintcf,nextci,nextcf and allocate needed
 * for computation arrays of needed size
 *
 * @param read_type
 * @param myrank
 * @param nprocs
 * @param nintci
 * @param nintcf
 * @param nextci
 * @param nextcf
 * @param lcc
 * @param points_count
 * @param points
 * @param elems
 * @param intcell_per_proc
 * @param extcell_per_proc
 * @param local_global_index
 * @return
 */
int allocate_local_variables(char* read_type, int myrank, int nprocs,
		int *nintci, int *nintcf, int *nextci,int *nextcf,
		int ***lcc,
		double **bs, double **be, double **bn, double **bw,
		double **bl, double **bh, double **bp, double **su,
		int* points_count, int*** points, int** elems,
		int **local_global_index,
		int *intcell_per_proc, int *extcell_per_proc,
		int *local_global_index_g, int points_count_g);

/**
 * Send and receive or read all needed local data
 *
 * @param read_type
 * @param myrank
 * @param nprocs
 * @param nintci
 * @param nintcf
 * @param nextci
 * @param nextcf
 * @param lcc
 * @param bs
 * @param be
 * @param bn
 * @param bw
 * @param bl
 * @param bh
 * @param bp
 * @param su
 * @param points_count
 * @param points
 * @param elems
 * @param intcell_per_proc
 * @param extcell_per_proc
 * @param local_global_index
 * @return
 */
int send_or_read_data(char* read_type, int myrank, int nprocs,
		int nintci, int nintcf, int nextci,int nextcf,
		int **lcc,
		double *bs, double *be, double *bn, double *bw,
		double *bl, double *bh, double *bp, double *su,
		int points_count, int** points, int* elems,
		int *local_global_index,
		int *intcell_per_proc, int *extcell_per_proc,
		int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
		int **lcc_g,
		double **bs_g, double **be_g, double **bn_g, double **bw_g,
		double **bl_g, double **bh_g, double **bp_g, double **su_g,
		int points_count_g, int** points_g, int **elems_g,
		int *local_global_index_g);

// TODO: write comment + rename ne here and in *.c
void fill_local_global_index(int nprocs,int *local_global_index_g, int ne,  int *metis_idx, int *intcell_per_proc);

// TODO: comment
int sort_data_by_local_global_index(int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
	double **bs_g, double **be_g, double **bn_g, double **bw_g,
	double **bl_g, double **bh_g, double **bp_g, double **su_g,
	int **elems_g,
	int *local_global_index_g);

/**
 * Counts the number of externall cell for each process and change lcc to have
 * indexes which are local
 *
 * TODO:better comment and rename the function
 */
void count_ext_cells(int nprocs, int *local_global_index_g,
		int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
		int **lcc_g, int *metis_idx,
		int *intcell_per_proc, int *extcell_per_proc);

#endif /* INITIALIZATION_ALGORITHMS_H_ */
