#ifndef FINALIZATION_H_
#define FINALIZATION_H_

void finalization(char* file_in, int nprocs, int myrank, int total_iters, double residual_ratio,
				int nintci, int nintcf, double* var, int* local_global_index, int* global_local_index);

#endif /* FINALIZATION_H_ */
