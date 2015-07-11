#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
	int i;
	char *file_in_name = argv[1], *file_out_name = argv[2];
    int nintci, nintcf;    /// internal cells start and end index
    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci, nextcf;
    int *lcc;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients
    double *boundaries;

    FILE *fp_in = fopen(file_in_name, "r");
    if (fp_in == NULL) {
        printf("Error opening file %s to read\n", file_in_name);
        return -1;
    }
    FILE *fp_out=fopen(file_out_name, "wb");
    if (fp_out == NULL) {
        printf("Error opening file %s to write\n", file_in_name);
        fclose(fp_in);
        return -1;
    }
    //4 variables in total!!!
    fscanf(fp_in, "%d", &nintci);
    fwrite(&nintci, sizeof(nintci), 1, fp_out);
    fscanf(fp_in, "%d", &nintcf);
    fwrite(&nintcf, sizeof(nintcf), 1, fp_out);
    fscanf(fp_in, "%d", &nextci);
    fwrite(&nextci, sizeof(nextci), 1, fp_out);
    fscanf(fp_in, "%d", &nextcf);
    fwrite(&nextcf, sizeof(nextcf), 1, fp_out);
    //allocating lcc
    if ((lcc = (int*) malloc( 6 * ( (nintcf) - (nintci) + 1) * sizeof(int) )) == NULL) {
        printf("malloc failed to allocate lcc");
        return -1;
    }
    // start reading lcc
    for (i = 0; i < 6*(nintcf-nintci+1); i++) {
		fscanf(fp_in, "%d", &(lcc)[i]);
    }
    // write whole array
	fwrite(lcc, sizeof(int), ( (nintcf) - (nintci) + 1) * 6, fp_out);
    // allocate other arrays
    if ((boundaries = (double *) malloc(8*(nintcf - nintci + 1) * sizeof(double))) == NULL) {
        printf("malloc() failed to allocate boundaries\n");
        return -1;
    }
    // read and write boundary information
    for (i = 0; i < 8*(nintcf - nintci + 1); i++) {
        fscanf(fp_in, "%lf", &((boundaries)[i]));
    }
    fwrite(boundaries, sizeof(double), 8*(nintcf - nintci + 1) , fp_out);
    // Finalize
    free(boundaries);
    free(lcc);
    fclose(fp_in);
    fclose(fp_out);
    printf("File '%s' converted successfully in '%s'\n",file_in_name, file_out_name);
	return 0;
}
