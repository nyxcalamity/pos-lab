/*
 * initialization_algorithms.h
 *
 *  Created on: Nov 10, 2014
 *      Author: power-morzh
 */

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
int compute_metis(char* part_type, char* read_type, int myrank, int nprocs, int nintci_g,
        int nintcf_g, int nextci_g, int nextcf_g, int *nintci, int *nintcf, int *nextci, int *nextcf,
        int **lcc_g, int points_count_g, int**points_g, int* elems_g, int *intcell_per_proc,
        int *extcell_per_proc, int** local_global_index_g, int** local_global_index, int **metis_idx);

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
int allocate_lcc_elems_points(char* read_type, int myrank, int nprocs, int *nintci, int *nintcf,
       int ***lcc, int* points_count, int*** points,
        int** elems, int **local_global_index, int points_count_g);

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
int fill_lcc_elems_points(char* read_type, int myrank, int nprocs, int nintci, int nintcf,
        int **lcc, int points_count, int** points, int* elems,
        int *local_global_index,  int **lcc_g,
        int points_count_g, int** points_g, int **elems_g);

int allocate_boundary_coef(char* read_type, int myrank, int nprocs, int *nextcf, double **bs, double **be, double **bn, double **bw,
        double **bl, double **bh, double **bp, double **su);

int fill_boundary_coef(char* read_type, int myrank, int nprocs, int nintci, int nintcf, int nextci,
        int nextcf, double *bs, double *be, double *bn, double *bw, double *bl,
        double *bh, double *bp, double *su,
        int *local_global_index, double **bs_g, double **be_g,
        double **bn_g, double **bw_g, double **bl_g, double **bh_g, double **bp_g, double **su_g);

void fill_local_global_index(char* read_type, int myrank, int nintci, int nintcf,
        int** local_global_index, int *metis_idx, int nelems_g);

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
 * FIXME: delete not needed arguments
 */
void build_lists_g2l_next(char* part_type, char* read_type, int nprocs, int myrank,
        int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst,
        int **recv_cnt, int*** recv_lst);

/**
 * Sends and receives sizes of send list(which is equal to size of receive list) and
 * allocate needed memory
 *
 * FIXME: delete unused arguments + better comments
 */
void allocate_send_lists(char* part_type, char* read_type, int nprocs, int myrank,
        int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst,
        int **recv_cnt, int*** recv_lst);

/**
 * Send receive lists and get them in send lists.
 *
 * FIXME: better comment + delete unused arguments
 */
void exchange_lists(char* part_type, char* read_type, int nprocs, int myrank,
        int* nintci, int* nintcf, int* nextci,
        int* nextcf, int*** lcc, int* points_count,
        int*** points, int** elems, double** var, double** cgup, double** oc,
        double** cnorm, int** local_global_index, int** global_local_index,
        int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst,
        int **recv_cnt, int*** recv_lst);

#endif /* INITIALIZATION_ALGORITHMS_H_ */
