#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_

int initialization(char* file_in, int input_key, int part_key, int read_key, int nprocs, int myrank, 
        int* nintci, int* nintcf, int* nextci, int* nextcf, int*** lcc, double** bs, double** be, 
        double** bn, double** bw, double** bl, double** bh, double** bp, double** su, int* points_count, 
        int*** points, int** elems, double** var, double** cgup, double** oc, double** cnorm, 
        int** local_global_index, int** global_local_index, int *nghb_cnt, int** nghb_to_rank, 
        int** send_cnt, int*** send_lst, int **recv_cnt, int*** recv_lst);

#endif /* INITIALIZATION_H_ */

