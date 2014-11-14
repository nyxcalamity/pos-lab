/*
 * initialization_algorithms.c
 *
 *  Created on: Nov 10, 2014
 *      Author: Denys Korzh, Denys Sobchyshak
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "initialization_algorithms.h"
#include "util_read_files.h"

#include "mpi.h"

#include "metis.h"

#define NINTCF_SEND_INDEX 0
#define NEXTCF_SEND_INDEX 1
#define POINTSCOUNT_SEND_INDEX 2

int read_global_data_or_geometry(char* file_in,char* read_type, int myrank,
        int *nintci, int *nintcf, int *nextci, int *nextcf,
        int ***lcc,
        double **bs, double **be, double **bn, double **bw,
		double **bl, double **bh, double **bp, double **su,
		int* points_count, int***points, int** elems) {
	int f_status=0;
	if ( !strcmp( read_type, "oneread" ) ) {
		if(myrank == 0) {
			// read-in the input file
			f_status = read_binary_geo(file_in, &*nintci, &*nintcf, &*nextci, &*nextcf, &*lcc, &*bs,
									   &*be, &*bn, &*bw, &*bl, &*bh, &*bp, &*su, &*points_count,
									   &*points, &*elems);
		}
	} else {
//TODO: Implement reading of geometry by each process
	}
    if ( f_status != 0 ) return f_status;
	return 0;
}


int compute_metis(char* part_type, char* read_type, int myrank, int nprocs,
        int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
		int *nintci, int *nintcf, int *nextci, int *nextcf,
        int **lcc,
		int points_count_g, int**points_g, int* elems_g,
		int *intcell_per_proc, int *extcell_per_proc, int** local_global_index_g,
		int **metis_idx) {
	int i=0;
	/** Metis variables **/
	// TODO: Think of better names
    idx_t ne, nn, ncommon, nparts, objval_idx;
    idx_t *eptr, *eind, *epart_idx, *npart_idx;
    int *npart;
    /** End Metis variables **/
	if( !strcmp( read_type, "oneread" ) && myrank==0 ) {
		// Allocate
		if ( (*local_global_index_g = (int *) malloc((nintcf_g + 1) * sizeof(int))) == NULL ) {
			fprintf(stderr, "malloc(local_global_index) failed\n");
			return -1;
		}

		if ( !strcmp( part_type, "classic" ) ) {

			for(i=0; i<nprocs-1; i++) {
				intcell_per_proc[i] = (nintcf_g + 1 + (nprocs-1))/nprocs;	// (nprocs-1) is important for not to loosing one cells in the processes
			}
			intcell_per_proc[nprocs-1] = nintcf_g + 1 -
					(nprocs-1)*((nintcf_g + 1 + (nprocs-1))/nprocs);  		//if our domain can't be divided in equal parts (breakets are very important!!!)
			for(i=0; i<nprocs-1; i++) {
				extcell_per_proc[i] = (nextcf_g - nextci_g + 1 + (nprocs-1))/nprocs;
			}
			extcell_per_proc[nprocs-1] = (nextcf_g - nextci_g + 1) -
					(nprocs-1)*((nextcf_g - nextci_g + 1 + (nprocs-1))/nprocs);  		//if our domain can't be divided in equal parts (breakets are very important!!!)
			for(i=0;i<(nintcf_g + 1);i++) {
				(*local_global_index_g)[i]=i;	// In classical it's simple
			}

		} else {
		    ne = nintcf_g - nintci_g + 1;
		    nn = points_count_g;
		    ncommon = 4;
		    nparts = nprocs;

		    eptr = ( idx_t* ) calloc( sizeof( idx_t ), ( ne + 1 ) );
		    eind = ( idx_t* ) calloc( sizeof( idx_t ), ( ne * 8 ) );
		    *metis_idx = ( int* ) calloc( sizeof( int ), ( ne ) );
		    npart = ( int* ) calloc( sizeof( int ), ( nn ) );
		    epart_idx = ( idx_t* ) calloc( sizeof( idx_t ), ( ne ) );
		    npart_idx = ( idx_t* ) calloc( sizeof( idx_t ), ( nn ) );

		    for( i = 0; i < ( ne + 1 ); i++ ) {
		        eptr[i] = 8 * i;
		    }
		    for( i = 0; i < ( ne * 8 ); i++ ) {
		         eind[i] = elems_g[i];
		    }
		    if ( !strcmp( part_type, "classic" ) ) {
				METIS_PartMeshDual( &ne,
						 &nn,
						 eptr,
						 eind,
						 NULL,
						 NULL,
						 &ncommon,
						 &nparts,
						 NULL,
						 NULL,
						 &objval_idx,
						 epart_idx,
						 npart_idx );
		    } else {
				METIS_PartMeshNodal( &ne,
						 &nn,
						 eptr,
						 eind,
						 NULL,
						 NULL,
						 &nparts,
						 NULL,
						 NULL,
						 &objval_idx,
						 epart_idx,
						 npart_idx );
			}
	        for(i = 0; i < ne; i++) {
	        	(*metis_idx)[i] = ( int )epart_idx[i];
	        }
		    for(i=0;i<20;i++) {
		    	printf("!!!!!!!%d\n", (*metis_idx)[i]);
		    }
			for(i=0; i<nprocs; i++) {
				intcell_per_proc[i] = 0;
			}
			for(i=0; i<nprocs-1; i++) {
				extcell_per_proc[i] = (nextcf_g - nextci_g + 1 + (nprocs-1))/nprocs;
			}
			extcell_per_proc[nprocs-1] = (nextcf_g - nextci_g + 1) -
					(nprocs-1)*((nextcf_g - nextci_g + 1 + (nprocs-1))/nprocs);  		//if our domain can't be divided in equal parts (breakets are very important!!!)
			fill_local_global_index(nprocs,*local_global_index_g, ne, *metis_idx,intcell_per_proc);
		}
	}
	return 0;
}

int allocate_local_variables(char* read_type, int myrank, int nprocs,
		int *nintci, int *nintcf, int *nextci,int *nextcf,
		int ***lcc,
		double **bs, double **be, double **bn, double **bw,
		double **bl, double **bh, double **bp, double **su,
		int* points_count, int*** points, int** elems,
		int **local_global_index,
		int *intcell_per_proc, int *extcell_per_proc,
		int *local_global_index_g, int points_count_g) {
	int i=0;
	MPI_Status status ;
	if( !strcmp( read_type, "oneread" ) ) {
		// Before we allocate, we need to know how much memory to allocate
		*nintci=0;	// it is the same for all processes
		if(myrank==0) {
			*nintcf = intcell_per_proc[0]-1;	// why to communicate if we already have it
			*nextci = intcell_per_proc[0];
			*nextcf = *nintcf + extcell_per_proc[0];
			*points_count = points_count_g;
			for(i=1; i<nprocs; i++) {
				MPI_Send(&intcell_per_proc[i],1,MPI_INT,i,NINTCF_SEND_INDEX,MPI_COMM_WORLD);
				MPI_Send(&extcell_per_proc[i],1,MPI_INT,i,NEXTCF_SEND_INDEX,MPI_COMM_WORLD);
				MPI_Send(&points_count_g,1,MPI_INT,i,POINTSCOUNT_SEND_INDEX,MPI_COMM_WORLD);
			}
		} else {
			MPI_Recv(nintcf,1, MPI_INT, 0, NINTCF_SEND_INDEX, MPI_COMM_WORLD, &status);
			MPI_Recv(nextcf,1, MPI_INT, 0, NEXTCF_SEND_INDEX, MPI_COMM_WORLD, &status);
			MPI_Recv(points_count,1, MPI_INT, 0, POINTSCOUNT_SEND_INDEX, MPI_COMM_WORLD, &status);
			*nextci = *nintcf;
			*nintcf = *nintcf-1;	// We need to subtract, because we got the number of cells(not the last index!)
			*nextcf = *nintcf + *nextcf;	// We need to subtract, because we got the number of cells(not the last index!)
		}
	}
	// Allocation
    if ( (*lcc = (int**) malloc((*nintcf + 1) * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }
    for ( i = 0; i < *nintcf + 1; i++ ) {
        if ( ((*lcc)[i] = (int *) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of lcc\n");
            return -1;
        }
    }
    if ( (*elems = (int *) malloc((*nintcf + 1) * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc(elems) failed\n");
        return -1;
    }
    if ( (*local_global_index = (int *) malloc((*nintcf + 1) * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc(local_global_index) failed\n");
        return -1;
    }
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
    if ( (*points = (int **) calloc(*points_count, sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }
    for ( i = 0; i < *points_count; i++ ) {
        if ( ((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }
    // End allocation
return 0;
}


int send_or_read_data(char* read_type, int myrank, int nprocs,
		int nintci, int nintcf, int nextci,int nextcf,
		int **lcc,
		double *bs, double *be, double *bn, double *bw,
		double *bl, double *bh, double *bp, double *su,
		int points_count, int** points, int* elems,
		int *local_global_index,
		int *intcell_per_proc, int *extcell_per_proc,
		int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
		int **lcc_g,
		double **bs_g, double **be_g, double **bn_g, double **bw_g,
		double **bl_g, double **bh_g, double **bp_g, double **su_g,
		int points_count_g, int** points_g, int **elems_g,
		int *local_global_index_g) {
	int k = 0;
	int i=0;
	MPI_Status status ;
	if( !strcmp( read_type, "oneread" ) ) {
		if(myrank==0) {

			// TODO: check on bad allocation
			sort_data_by_local_global_index(nintci_g, nintcf_g, nextci_g, nextcf_g,
					lcc_g,
					&*bs_g, &*be_g, &*bn_g, &*bw_g,
					&*bl_g, &*bh_g, &*bp_g, &*su_g,
					&*elems_g,
					local_global_index_g);

			int start_idx = 0;	// Used to send data, from some index(because nodes can have not equal number of cells)
			// Copy memory for process 0
			for(i=nintci; i<nintcf + 1; i++) {
				memcpy(lcc[i], lcc_g[i], 6*sizeof(int));
			}
			memcpy(bs, *bs_g, intcell_per_proc[0] * sizeof(double));
			memcpy(be, *be_g, intcell_per_proc[0] * sizeof(double));
			memcpy(bn, *bn_g, intcell_per_proc[0] * sizeof(double));
			memcpy(bw, *bw_g, intcell_per_proc[0] * sizeof(double));
			memcpy(bl, *bl_g, intcell_per_proc[0] * sizeof(double));
			memcpy(bh, *bh_g, intcell_per_proc[0] * sizeof(double));
			memcpy(bp, *bp_g, intcell_per_proc[0] * sizeof(double));
			memcpy(su, *su_g, intcell_per_proc[0] * sizeof(double));
			for ( i = 0; i < points_count; i++ ) {
				memcpy(points[i], points_g[i], 3*sizeof(int));
			}
			memcpy(elems, *elems_g, intcell_per_proc[0]*8*sizeof(int));
			memcpy(local_global_index, local_global_index_g,
					intcell_per_proc[0] * sizeof(int));
			// Send all other data
			for (k=1; k<nprocs; ++k) {
				start_idx += intcell_per_proc[k-1];
				for(i=0; i < intcell_per_proc[k]; ++i) {
					MPI_Send(lcc_g[start_idx + i],6,MPI_INT,k,
							0,MPI_COMM_WORLD);
				}
				MPI_Send(&(*bs_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 1, MPI_COMM_WORLD);
				MPI_Send(&(*be_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 2, MPI_COMM_WORLD);
				MPI_Send(&(*bn_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 3, MPI_COMM_WORLD);
				MPI_Send(&(*bw_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 4, MPI_COMM_WORLD);
				MPI_Send(&(*bl_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 5, MPI_COMM_WORLD);
				MPI_Send(&(*bh_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 6, MPI_COMM_WORLD);
				MPI_Send(&(*bp_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 7, MPI_COMM_WORLD);
				MPI_Send(&(*su_g)[start_idx],intcell_per_proc[k],
						MPI_DOUBLE, k, 8, MPI_COMM_WORLD);
				for ( i = 0; i < points_count; i++ ) {
					MPI_Send(points_g[i],3,MPI_INT,k,99,MPI_COMM_WORLD);
				}
				MPI_Send(&(*elems_g)[start_idx*8],intcell_per_proc[k]*8,
						MPI_INT, k, 100, MPI_COMM_WORLD);
				MPI_Send(&local_global_index_g[start_idx],
						intcell_per_proc[k], MPI_INT, k, 101, MPI_COMM_WORLD);

			}
		} else {
			// Receive data from process 0
			for(i=nintci; i < nintcf; ++i) {
				MPI_Recv(lcc[i],6, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			}
			MPI_Recv(bs,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(be,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
			MPI_Recv(bn,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
			MPI_Recv(bw,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);
			MPI_Recv(bl,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status);
			MPI_Recv(bh,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &status);
			MPI_Recv(bp,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, &status);
			MPI_Recv(su,(nintcf-nintci+1),
					MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, &status);
			for ( i = 0; i < points_count; i++ ) {
				MPI_Recv(points[i],3, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);
			}
			MPI_Recv(elems,(nintcf-nintci+1)*8,
					MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(local_global_index,(nintcf-nintci+1),
					MPI_INT, 0, 101, MPI_COMM_WORLD, &status);
		}
	}
	return 0;
}

// TODO: change name of ne
void fill_local_global_index(int nprocs,int *local_global_index_g, int ne,  int *metis_idx, int *intcell_per_proc) {
	int i=0;
	int k=0;
	int writer_counter=0;
	for (k=0;k<nprocs;++k) {
		for(i=0; i<ne; ++i) {
			if(metis_idx[i]==k) {
				++intcell_per_proc[k];
				local_global_index_g[writer_counter]=i;
				++writer_counter;
			}
		}
	}
}

int sort_data_by_local_global_index(int nintci_g, int nintcf_g, int nextci_g, int nextcf_g,
		int **lcc_g,
		double **bs_g, double **be_g, double **bn_g, double **bw_g,
		double **bl_g, double **bh_g, double **bp_g, double **su_g,
		int **elems_g,
		int *local_global_index_g) {
	int i=0;
	int *tmp_elems;	// This two arrays are needed to send data in correct order
	double *tmp_b;
	void *tmp; // Not to lose the pointer
	if ((tmp_elems = (int *) malloc((nintcf_g+1)*8 * sizeof(int))) == NULL)
	{
		printf("malloc() of tmp_elems in send_or_read_data failed\n");
		return -1;
	}
	if ((tmp_b = (double *) malloc((nintcf_g + 1) * sizeof(double))) == NULL)
	{
		printf("malloc() of tmp_b in send_or_read_data failed\n");
		return -1;
	}
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*bs_g)[local_global_index_g[i]];
		memcpy(&tmp_elems[8*i], &(*elems_g)[8*local_global_index_g[i]],8*sizeof(int));
	}
//	memcpy(elems_g,tmp_elems,8*(nintcf_g+1)*sizeof(int));
	tmp = *elems_g;
	*elems_g = tmp_elems;
	tmp_elems = tmp;
	tmp=*bs_g;
	*bs_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*be_g)[local_global_index_g[i]];
	}
	tmp=*be_g;
	*be_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*bn_g)[local_global_index_g[i]];
	}
	tmp=*bn_g;
	*bn_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*bw_g)[local_global_index_g[i]];
	}
	tmp=*bw_g;
	*bw_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*bl_g)[local_global_index_g[i]];
	}
	tmp=*bl_g;
	*bl_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*bh_g)[local_global_index_g[i]];
	}
	tmp=*bh_g;
	*bh_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*bp_g)[local_global_index_g[i]];
	}
	tmp=*bp_g;
	*bp_g=tmp_b;
	tmp_b=tmp;
	for(i=0; i<nintcf_g+1; ++i) {
		tmp_b[i] = (*su_g)[local_global_index_g[i]];
	}
	tmp=*su_g;
	*su_g=tmp_b;
	tmp_b=tmp;
	free(tmp_b);
	free(tmp_elems);
	return 0;
}
