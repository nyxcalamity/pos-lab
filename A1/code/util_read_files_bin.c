#include <stdio.h>
#include <stdlib.h>
#include "util_read_files_bin.h"

int read_formatted_bin(char *filename, int *nintci, int *nintcf, int *nextci,
        int *nextcf, int ***lcc, double **bs, double **be, double **bn,
        double **bw, double **bl, double **bh, double **bp, double **su)
{
    int i;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("Error opening file %s\n", filename);
        return -1;
    }
    //4 variables in total!!!
    fread(nintci, sizeof(*nintci),1,fp);
    fread(nintcf, sizeof(*nintcf),1,fp);
    fread(nextci, sizeof(*nextci),1,fp);
    fread(nextcf, sizeof(*nextcf),1,fp);

    //allocating lcc
    if ((*lcc = (int**) malloc( ( (*nintcf) - (*nintci) + 1) * sizeof(int*) )) == NULL)
    {
        printf("malloc failed to allocate first dimension of lcc (nintcf)");
        return -1;
    }
    for (i = (*nintci); i <= (*nintcf); i++)
    {
        if (((*lcc)[i] = (int *) malloc( 6 * sizeof(int) )) == NULL)
        {
            printf("malloc(lcc) failed\n");
            return -1;
        }
    }

    //start reading lcc
    //note that c array index starts from 0 while fortran starts from 1!
    for (i = (*nintci); i <= (*nintcf); i++)
    {
    	fread((*lcc)[i], sizeof(int),6,fp);
    }

    // allocate other arrays
    if ((*bs = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*be = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bn = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bw = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bl = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bh = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bp = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*su = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }

    // read the other arrays
    for (i = (*nintci); i <= *nintcf; i++)
    {
        fread(&(*bs)[i], sizeof((*bs)[i]), 1, fp);
        fread(&(*be)[i], sizeof((*be)[i]), 1, fp);
        fread(&(*bn)[i], sizeof((*bn)[i]), 1, fp);
        fread(&(*bw)[i], sizeof((*bw)[i]), 1, fp);
        fread(&(*bl)[i], sizeof((*bl)[i]), 1, fp);
        fread(&(*bh)[i], sizeof((*bh)[i]), 1, fp);
        fread(&(*bp)[i], sizeof((*bp)[i]), 1, fp);
        fread(&(*su)[i], sizeof((*su)[i]), 1, fp);
    }

    fclose(fp);

    return 0;
}
