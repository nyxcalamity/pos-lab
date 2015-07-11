#ifndef COMPUTE_SOLUTION_H_
#define COMPUTE_SOLUTION_H_

int compute_solution(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int nghb_cnt,
                     int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst,
                     char *file_in, int points_count, int **points, int *elems, char* part_type, char* read_type,
                     int **l2g_g, int *int_cells_per_proc);

#endif /* COMPUTE_SOLUTION_H_ */

