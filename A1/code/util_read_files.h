#ifndef XREAD_H_
#define XREAD_H_

int read_formatted(char *filename,
                int *nintci, int *nintcf, int *nextci, int *nextcf,
                int ***lcc,
                double **bs, double **be, double **bn, double **bw, double **bl, double **bh, double **bp, double **su);
#endif /* XREAD_H_ */


