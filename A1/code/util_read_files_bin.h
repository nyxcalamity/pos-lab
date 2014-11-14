/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#ifndef XREADBIN_H_
#define XREADBIN_H_

int read_formatted_bin(char *filename,
                int *nintci, int *nintcf, int *nextci, int *nextcf,
                int ***lcc,
                double **bs, double **be, double **bn, double **bw, double **bl, double **bh, double **bp, double **su);
#endif /* XREADBIN_H_ */


