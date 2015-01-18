/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util_write_files.h"

int compute_solution(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int nghb_cnt,
                     int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst,
                     char *file_in, int points_count, int **points, int *elems, char* part_type, char* read_type,
                     int **l2g_g, int *int_cells_per_proc) {
    int rank=0;
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

    /** the reference residual */
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }

    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }

    /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
//            if(nc==36507)
//                printf("direc1=%.15lf,resvec=%.15lf,cgup=%.15lf\n",direc1[nc],resvec[nc],cgup[nc]);
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
//            if(nc==36507)
//                printf("direc2=%.15lf,resvec=%.15lf,cgup=%.15lf\n",direc2[nc],resvec[nc],cgup[nc]);
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
                         - be[nc] * direc1[lcc[nc][1]] - bn[nc] * direc1[lcc[nc][2]]
                         - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                         - bh[nc] * direc1[lcc[nc][5]];
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        double occ_a[nprocs];
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

//            for ( nc = nintci; nc <= nintcf; nc++ ) {
//                occ = occ + direc2[nc] * adxor1[nc];
//            }
            // This sum we did to imitate sum with MPI_Allreduce and we got the same result.
            for (rank=0; rank<nprocs; ++rank) {
                occ_a[rank]=0;
                for ( nc = 0; nc < int_cells_per_proc[rank]; nc++ ) {
                    occ_a[rank] = occ_a[rank] + direc2[l2g_g[rank][nc]] * adxor1[l2g_g[rank][nc]];
                }
                occ = occ + occ_a[rank];
            }

            oc1 = occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

//                for ( nc = nintci; nc <= nintcf; nc++ ) {
//                    occ = occ + direc2[nc] * adxor1[nc];
//                }
                // This sum we did to imitate sum with MPI_Allreduce and we got the same result.
                for (rank=0; rank<nprocs; ++rank) {
                    occ_a[rank]=0;
                    for ( nc = 0; nc < int_cells_per_proc[rank]; nc++ ) {
                        occ_a[rank] = occ_a[rank] + direc2[l2g_g[rank][nc]] * adxor1[l2g_g[rank][nc]];
                    }
                    occ = occ + occ_a[rank];
                }

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
//                for ( nc = nintci; nc <= nintcf; nc++ ) {
//                    occ = occ + direc2[nc] * adxor2[nc];
//                }
                // This sum we did to imitate sum with MPI_Allreduce and we got the same result.
                for (rank=0; rank<nprocs; ++rank) {
                    occ_a[rank]=0;
                    for ( nc = 0; nc < int_cells_per_proc[rank]; nc++ ) {
                        occ_a[rank] = occ_a[rank] + direc2[l2g_g[rank][nc]] * adxor2[l2g_g[rank][nc]];
                    }
                    occ = occ + occ_a[rank];
                }

                oc2 = occ / cnorm[2];
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }
//        printf("i%d, occ=%.15lf\n",iter,occ);

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
//        for ( nc = nintci; nc <= nintcf; nc++ ) {
//            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
//            omega = omega + resvec[nc] * direc2[nc];
//        }
        // This sum we did to imitate sum with MPI_Allreduce and we got the same result.
        double omega_a[nprocs];
        double cnorm_a[nprocs];
        for (rank=0; rank<nprocs; ++rank) {
            omega_a[rank]=0;
            cnorm_a[rank]=0;
            for ( nc = 0; nc < int_cells_per_proc[rank]; nc++ ) {
                cnorm_a[rank] = cnorm_a[rank] + direc2[l2g_g[rank][nc]] * direc2[l2g_g[rank][nc]];
                omega_a[rank] = omega_a[rank] + resvec[l2g_g[rank][nc]] * direc2[l2g_g[rank][nc]];
            }
            omega = omega + omega_a[rank];
            cnorm[nor] = cnorm[nor] + cnorm_a[rank];
//            printf("r%d, omega=%.30lf\n",rank, omega_a[rank]);
//            printf("r%d, cnorm[nor]=%.30lf\n",rank, cnorm_a[rank]);
        }
//        printf("r%d, cnorm[nor]=%.30lf\n",myrank, cnorm[nor]);
//        printf("r%d, omega=%.30lf\n",myrank, omega);

        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
//            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }
//        printf("r%d, omega=%.20lf\n",myrank, omega);
        // This sum we did to imitate sum with MPI_Allreduce and we got the same result.
        double res_updated_a[nprocs];
        for (rank=0; rank<nprocs; ++rank) {
            res_updated_a[rank]=0;
            for ( nc = 0; nc < int_cells_per_proc[rank]; nc++ ) {
                res_updated_a[rank] = res_updated_a[rank] + resvec[l2g_g[rank][nc]] * resvec[l2g_g[rank][nc]];
            }
            res_updated = res_updated + res_updated_a[rank];
        }
        res_updated = sqrt(res_updated);
        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
        if (iter==max_iters) {
            check_compute_values(file_in, part_type, read_type, nprocs, myrank,
                    nintci, nintcf, nextcf, omega, nor,
                    resvec, direc1, direc2, var, cnorm, l2g_g, int_cells_per_proc);
        }
    }
//    vtk_check(file_in, myrank, nintci, nintcf, resvec, direc1, direc2, var, points_count, points,
//                    elems, local_global_index, (nintcf-nintci+1));

    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);
    free(resvec);

    return iter;
}


