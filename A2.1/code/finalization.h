#ifndef FINALIZATION_H_
#define FINALIZATION_H_

void finalization(char* file_in, char* out_prefix, int total_iters, double residual_ratio,
                  int nintci, int nintcf, int points_count, int** points, int* elems, double* var,
                  double* cgup, double* su);

#endif /* FINALIZATION_H_ */

