#ifndef XREAD_H_
#define XREAD_H_

int read_binary_geo(char *file_name, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
                    double **BS, double **BE, double **BN, double **BW, double **BL, double **BH,
                    double **BP, double **SU, int* nodeCnt, int***points, int** elems);


int read_lcc_boundary(char *file_name, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
                    double **BS, double **BE, double **BN, double **BW, double **BL, double **BH,
                    double **BP, double **SU);


int read_geometry(char *file_name, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF,
        int* points_count, int*** points, int** elems);


int read_lcc_local(char *file_name, int NINTCI, int NINTCF, int **LCC, int *local_global_index);


int read_boundary_local(char *file_name, int NINTCI, int NINTCF, int nintcf_g,
                    double *BS, double *BE, double *BN, double *BW, double *BL, double *BH,
                    double *BP, double *SU, int *local_global_index);

#endif /* XREAD_H_ */