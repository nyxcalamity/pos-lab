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
#include <mpi.h>


int compute_solution(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int nghb_cnt,
                     int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst) {
//    check_compute_arguments(nprocs, myrank, max_iters, nintci, nintcf, nextcf,
//                        lcc, bp, bs, bw, bl, bn, be, bh,
//                         cnorm, var, su, cgup, residual_ratio,
//                         local_global_index, global_local_index, nghb_cnt,
//                         nghb_to_rank, send_cnt, send_lst, recv_cnt, recv_lst,
//                         file_in, points_count, points, elems, part_type, read_type);
    MPI_Status status;
    MPI_Request request_send[nghb_cnt], request_recv[nghb_cnt];
    MPI_Datatype index_type[nghb_cnt];
    int *mpi_block_length[nghb_cnt], *mpi_displacements[nghb_cnt];
    
    /** buffers used to resend direc1 */
    int nghb_idx=0;

    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

    /** the reference residual */
    double resref_g = 0.0;
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }

    // Exchange resref
    MPI_Allreduce(&resref, &resref_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    resref=resref_g;
    resref = sqrt(resref);
    if (myrank == 0) {
        if ( resref < 1.0e-15 ) {
            fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
            MPI_Abort(MPI_COMM_WORLD, myrank);
            return 0;
        }
    }

    /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    
    //initialize mpi indexed data type arrays
    for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
        if ((mpi_block_length[nghb_idx] = (int *) malloc(send_cnt[nghb_idx]*sizeof(int))) == NULL) {
            fprintf(stderr, "malloc(mpi_block_length) failed\n");
            MPI_Abort(MPI_COMM_WORLD, myrank);
            return -1;
        }
        if ((mpi_displacements[nghb_idx] = (int *) malloc(send_cnt[nghb_idx]*sizeof(int))) == NULL) {
            fprintf(stderr, "malloc(mpi_displacements) failed\n");
            MPI_Abort(MPI_COMM_WORLD, myrank);
            return -1;
        }
    }

    while ( iter < max_iters ) {
//        if (iter==max_iters-1) {
//            check_compute_values(file_in, part_type, read_type, nprocs, myrank,
//                    nintci, nintcf, nextcf,
//                    resvec, direc1, direc2, var, cnorm);
//        }
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
//            if(local_global_index[nc]==36507)
//                printf("direc1=%.15lf,resvec=%.15lf,cgup=%.15lf\n",direc1[nc],resvec[nc],cgup[nc]);
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }       

        /** START Exchange of direc1 with MPI **/
        for(nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
            // Receive direc1
            MPI_Irecv(&direc1[recv_lst[nghb_idx][0]], recv_cnt[nghb_idx], MPI_DOUBLE, nghb_to_rank[nghb_idx],
                    nghb_to_rank[nghb_idx], MPI_COMM_WORLD, &request_recv[nghb_idx]);

            //initialize and register mpi data type
            for(nc=0; nc<send_cnt[nghb_idx]; ++nc) {
                mpi_block_length[nghb_idx][nc] = 1;
                mpi_displacements[nghb_idx][nc] = send_lst[nghb_idx][nc];
            }
            MPI_Type_indexed(send_cnt[nghb_idx], mpi_block_length[nghb_idx], mpi_displacements[nghb_idx], 
                MPI_DOUBLE, &index_type[nghb_idx]);
            MPI_Type_commit(&index_type[nghb_idx]);
            
            // Send direc1
            MPI_Isend(direc1, 1, index_type[nghb_idx], nghb_to_rank[nghb_idx],
                    myrank, MPI_COMM_WORLD, &request_send[nghb_idx]);
        }
        // Synchronize everything
        for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
            MPI_Wait(&request_recv[nghb_idx], &status);
            MPI_Wait(&request_send[nghb_idx], &status);
        }
//        if (iter==max_iters-1) {
//            check_compute_values(file_in, part_type, read_type, nprocs, myrank,
//                    nintci, nintcf, nextcf,
//                    resvec, direc1, direc2, var, cnorm);
//        }
        /** STOP Exchange of direc1 with MPI **/
        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
                         - be[nc] * direc1[lcc[nc][1]] - bn[nc] * direc1[lcc[nc][2]]
                         - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                         - bh[nc] * direc1[lcc[nc][5]];
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        double occ_g;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;
            occ_g=0;

            for ( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }
            MPI_Allreduce(&occ, &occ_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            occ = occ_g;

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
                occ_g=0;

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor1[nc];
                }
                MPI_Allreduce(&occ, &occ_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                occ = occ_g;

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
                occ_g=0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor2[nc];
                }
                MPI_Allreduce(&occ, &occ_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                occ = occ_g;

                oc2 = occ / cnorm[2];
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        double omega_g = 0;
        double cnorm_g = 0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }
        MPI_Allreduce(&cnorm[nor], &cnorm_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        cnorm[nor] = cnorm_g;

        MPI_Allreduce(&omega, &omega_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        printf("r%d, omega_g=%.30lf\n",myrank, omega_g);
//        printf("r%d, cnorm[nor]=%.30lf\n",myrank, cnorm[nor]);
        omega = omega_g;
        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        double res_updated_g = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }
        MPI_Allreduce(&res_updated, &res_updated_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        res_updated = res_updated_g;
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
    }

    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);
    free(resvec);
    
    printf("[INFO] Completed compute_solution on task #%d\n", myrank);
    
    return iter;
}