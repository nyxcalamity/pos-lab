/**
 * Helper functions for writing results to VTK and text files
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util_write_files.h"
#include "test_functions.h"
#include "posl_definitions.h"

int store_simulation_stats(char *in_file_name, char *out_file_name, int nintci, int nintcf,
                           double *var, int total_iters, double residual_ratio) {
    double *points = (double *) malloc((nintcf + 1) * sizeof(double));
    int i1, i2, i3, i4, i5;

    int counter;
    for ( counter = nintci; counter <= nintcf; counter++ )
        points[counter] = counter;

    if ( nintcf <= 1 ) {
        fprintf(stderr, "Error: NINTCF <= 1\n");
        return -1;
    }

    i1 = nintcf + 1;

    while ( i1 != 0 ) {
        i1 = i1 / 2;
        i2 = nintcf + 1 - i1;
        i4 = 1;

        do {
            i3 = i4;
            do {
                i5 = i3 + i1;
                if ( var[i3] <= var[i5] ) break;

                double z_dum = var[i3], i_dum = points[i3];

                var[i3] = var[i5];
                points[i3] = points[i5];
                var[i5] = z_dum;
                points[i5] = i_dum;
                i3 = i3 - i1;
            } while ( i3 >= 1 );
            i4++;
        } while ( i4 < i2 );
    }

    FILE *fp = fopen(out_file_name, "w");
    if ( fp == NULL ) {
        fprintf(stderr, "Error opening file %s for writing\n", out_file_name);
        return -1;
    }

    printf("========================================\n");
    printf("= AVL -  Linear Equation Solver - GCCG =\n");
    printf("========================================\n\n");
    printf("Input File:  %s\n", in_file_name);
    printf("Output File:  %s\n", out_file_name);
    printf("No. of Active Cells:  %d\n", nintcf);
    printf("Iterations Count: %d\n", total_iters);
    printf("Residual Ratio: %e\n", residual_ratio);
    printf("========================================\n\n");

    fprintf(fp, "========================================\n");
    fprintf(fp, "= AVL -  Linear Equation Solver - GCCG =\n");
    fprintf(fp, "========================================\n");
    fprintf(fp, "\n\n");
    fprintf(fp, "Input File:  %s\n", in_file_name);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Output File:  %s\n", out_file_name);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "No. of Active Cells:  %d\n", nintcf + 1);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Iterations Count: %d\n", total_iters);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Residual Ratio: %e\n", residual_ratio);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Addresses Solution (Minima)  \t Addresses Solution (Maxima)\n");
    fprintf(fp, "===========================         ===========================\n");

    int N;
    for ( N = 1; N <= 10; N++ )
        fprintf(fp, "%8.0lf \t %lf \t\t %8.0lf \t %lf\n", points[N], var[N], points[nintcf - N + 1],
                var[nintcf - N + 1]);

    fprintf(fp, "========================================\n");

    fclose(fp);

    free(points);
    return 0;
}

void vtk_write_unstr_grid_header(const char *experiment_name, const char *out_file_name,
                                 int start_index, int end_index, int points_count, int **points,
                                 int *elems) {
    int i, j;
    FILE *fp = NULL;
    // Total number of elements (cells)
    int elem_count = end_index - start_index + 1;

    fp = fopen(out_file_name, "w");
    if ( fp == NULL ) {
        fprintf(stderr, "Failed to open %s", out_file_name);
        return;
    }

    fprintf(fp, "# vtk DataFile Version 3.1\n");
    fprintf(fp, "%s\n", experiment_name);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    /*
     * The first line notes how many points there will be and in which format they'll be supplied (here "DOUBLE").
     * Each following line contains the xyz-coordinates of a point.
     * Based on this list, each point is assigned an id. The first point has id 0, the second point id 1, and so forth.
     */
    fprintf(fp, "POINTS %d DOUBLE\n", points_count);
    for ( i = 0; i < points_count; i++ )
        fprintf(fp, "%d %d %d\n", points[i][0], points[i][1], points[i][2]);

    fprintf(fp, "\n");

    /*
     * The first line notes how many cells there will be and how many numbers total will be supplied in the CELLS-block.
     * Each following cell line starts with a number saying how many point IDs are to be read in that line (here 8) followed by the list of those point IDs.
     * Based on this list, each cell is assigned an id. The first cell has id 0, the second id 1, etc.
     */
    fprintf(fp, "CELLS %d %d\n", elem_count, elem_count * 9);
    for ( i = 0; i < elem_count; i++ ) {
        fprintf(fp, "8 ");
        for ( j = 0; j < 8; j++ )
            fprintf(fp, "%d ", elems[8 * i + j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    /*
     * The first line gives how many cell-types are to be set (= the number of cells given in CELLS).
     * The following line is a list of cell-types that are assigned to all cells. The first cell is of type "11", so is the second cell and so forth.
     * Cell type "11" is VTK_VOXEL (volume cell) and requires 8 point coordinates to define the  of each cell (given in CELLS).
     */
    fprintf(fp, "CELL_TYPES %d\n", elem_count);
    for ( i = 0; i < elem_count; i++ )
        fprintf(fp, "11 ");
    fprintf(fp, "\n\n");

    /*
     * Denotes the beginning of the datasets that will be assigned to the cell.
     * POINT_DATA can be used to assign the datasets to the points instead of the cells.
     */
    fprintf(fp, "CELL_DATA %d\n", elem_count);
    fprintf(fp, "\n");

    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
}

void vtk_append_double(const char *out_file_name, const char *var_name, int start_index,
                       int end_index, double *values) {
    int i;
    FILE *fp = NULL;

    if ( (fp = fopen(out_file_name, "a")) == NULL ) {
        fprintf(stderr, "Failed to open %s", out_file_name);
        return;
    }

    /*
     * The first line gives the name of the dataset (variable) and its type (here "DOUBLE")
     * The second line selects the color table to use and it is usually "LOOKUP_TABLE default".
     * The following lines contain a value of the dataset per line
     */
    fprintf(fp, "SCALARS %s DOUBLE\n", var_name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( i = start_index; i <= end_index; i++ )
        fprintf(fp, "%f\n", values[i]);

    fprintf(fp, "\n");

    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
}

void vtk_append_integer(const char *out_file_name, const char *var_name, int start_index,
                        int end_index, int *values) {
    int i;
    FILE *fp = NULL;

    fp = fopen(out_file_name, "a");
    if ( fp == NULL ) {
        fprintf(stderr, "Failed to open %s", out_file_name);
        return;
    }

    /*
     * The first line gives the name of the dataset (variable) and its type (here "INT")
     * The second line selects the color table to use and it is usually "LOOKUP_TABLE default".
     * The following lines contain a value of the dataset per line
     */
    fprintf(fp, "SCALARS %s INT\n", var_name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( i = start_index; i <= end_index; i++ )
        fprintf(fp, "%d\n", values[i]);

    fprintf(fp, "\n");

    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
}


// start_of_student_code---------------------------------------------------------------------------------
void vtk_for_process(const char *file_in, const char *file_vtk_out, int nintci, int nintcf, 
        int points_count, int **points,  int *elems, int *local_global_index, int local_num_elems, 
        double *scalars) {
    // write vtk file
    vtk_write_unstr_grid_header(file_in, file_vtk_out, nintci, nintcf, points_count, points, elems);
    vtk_append_double(file_vtk_out, "SCALARS", nintci, nintcf, scalars);
    printf("Data VTK file succesfully generated! \n");
}


//TODO:include all execution setup identifiers in the file name
int vtk_check(char *file_in, int myrank, int nintci, int nintcf, double *su, double *cgup, 
        int points_count, int **points,  int *elems, int *local_global_index, int local_num_elems) {
    char szFileName[80];
    int i=0;
    double *scalars;
    //TODO:externalize this string
    const char *kOutputDirectoryName = "./out/";
    
    //find base file name
    char *data_file = strrchr(file_in,'/')+1;
    //strip data file base name
    data_file = strndup(data_file, strchr(data_file, '.')-data_file);
    
    if ((scalars = (double *) malloc((nintcf+1)*sizeof(double))) == NULL) {
        fprintf(stderr, "malloc(local_global_index) failed\n");
        return -1;
    }
    for (i=0; i<local_num_elems; i++){
        scalars[i] = cgup[i];
    }
    // Output CGUP
    sprintf(szFileName, "%s%s.cgup.rank%i.vtk", kOutputDirectoryName, data_file, myrank);
    test_distribution(file_in, szFileName, local_global_index, local_num_elems, scalars);
    sprintf(szFileName, "%s%s.cgup.cutted.rank%i.vtk", kOutputDirectoryName,data_file, myrank);
    vtk_for_process(file_in, szFileName, 0, local_num_elems-1, points_count, points, elems, 
            local_global_index, local_num_elems, scalars);
    // Output SU
    for (i=0; i<local_num_elems; i++){
        scalars[i] = su[i];
    }
    sprintf(szFileName, "%s%s.su.rank%i.vtk", kOutputDirectoryName, data_file, myrank);
    test_distribution(file_in, szFileName, local_global_index, local_num_elems, scalars);
    sprintf(szFileName, "%s%s.su.cutted.rank%i.vtk", kOutputDirectoryName, data_file, myrank);
    vtk_for_process(file_in, szFileName, 0, local_num_elems-1, points_count, points, elems, 
            local_global_index, local_num_elems, scalars);
    free(data_file);
    return 0;
}
// FIXME: add an argument which says what we want to show(recv,send,both)
void vtk_check_lists(char *file_in, int myrank,
        int *local_global_index, int local_num_elems,
        int nghb_cnt, int* nghb_to_rank, int* send_cnt, int** send_lst,
        int *recv_cnt, int** recv_lst, int output_style) {
    int start_from=0;
    int i=0, nghb_idx=0;
    int local_num_elems_big=local_num_elems;
    // Compute full size of all elements in processor
    for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
        if(output_style == VTK_ALL || output_style == VTK_RECV_LST) {
            local_num_elems_big += recv_cnt[nghb_idx];
        }
        if(output_style == VTK_ALL || output_style == VTK_SEND_LST) {
            local_num_elems_big += send_cnt[nghb_idx];
        }
    }

    int local_global_index_big[local_num_elems_big];
    double scalars[local_num_elems_big];
    // Fill local_global_index_big
    memcpy(local_global_index_big, local_global_index, local_num_elems*sizeof(int));
    start_from = local_num_elems;
    if(output_style == VTK_ALL || output_style == VTK_RECV_LST) {
        for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
            for (i=0; i<recv_cnt[nghb_idx]; ++i) {
                local_global_index_big[start_from + i] = recv_lst[nghb_idx][i];
            }
            start_from += recv_cnt[nghb_idx];
        }
    }
    if(output_style == VTK_ALL || output_style == VTK_SEND_LST) {
        for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
            for (i=0; i<send_cnt[nghb_idx]; ++i) {
                local_global_index_big[start_from + i] = send_lst[nghb_idx][i];
            }
            start_from += send_cnt[nghb_idx];
        }
    }

    // Fill scalars
    // Cells which belongs to process 0 will have value 1 to distinguish it from other cells
    // and other processors are incremented by one!
    for(i=0; i<local_num_elems; ++i) {
        scalars[i] = myrank+1;
    }
    start_from = local_num_elems;
    if(output_style == VTK_ALL || output_style == VTK_RECV_LST) {
        for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
            for (i=0; i<recv_cnt[nghb_idx]; ++i) {
                scalars[start_from + i] = nghb_to_rank[nghb_idx] + 1;
            }
            start_from += recv_cnt[nghb_idx];
        }
    }
    if(output_style == VTK_ALL || output_style == VTK_SEND_LST) {
        for (nghb_idx=0; nghb_idx<nghb_cnt; ++nghb_idx) {
            for (i=0; i<send_cnt[nghb_idx]; ++i) {
                scalars[start_from + i] = nghb_to_rank[nghb_idx] + 0.5;
            }
            start_from += send_cnt[nghb_idx];
        }
    }

    // Write VTK
    char szFileName[80];
    //TODO:externalize this string
    const char *kOutputDirectoryName = "./out/";

    //find base file name
    char *data_file = strrchr(file_in,'/')+1;
    //strip data file base name
    data_file = strndup(data_file, strchr(data_file, '.')-data_file);

    if(output_style == VTK_SEND_LST) {
        sprintf(szFileName, "%s%s.Send.rank%i.vtk", kOutputDirectoryName, data_file, myrank);
    } else if(output_style == VTK_RECV_LST) {
        sprintf(szFileName, "%s%s.Recv.rank%i.vtk", kOutputDirectoryName, data_file, myrank);
    } else {
        sprintf(szFileName, "%s%s.RecvandSend.rank%i.vtk", kOutputDirectoryName, data_file, myrank);
    }


    test_distribution(file_in, szFileName, local_global_index_big, local_num_elems_big, scalars);
}

void vtk_check_neighbour(char *file_in, int myrank,
        int *local_global_index, int local_num_elems,
        int nghb_cnt, int* nghb_to_rank, int* send_cnt, int** send_lst,
        int *recv_cnt, int** recv_lst, int output_style, int neighbour) {
    int start_from=0;
    int i=0;
    int local_num_elems_big=local_num_elems;
    // Compute full size of all elements in processor
    if(output_style == VTK_ALL || output_style == VTK_RECV_LST) {
        local_num_elems_big += recv_cnt[neighbour];
    }
    if(output_style == VTK_ALL || output_style == VTK_SEND_LST) {
        local_num_elems_big += send_cnt[neighbour];
    }

    int local_global_index_big[local_num_elems_big];
    double scalars[local_num_elems_big];
    // Fill local_global_index_big
    memcpy(local_global_index_big, local_global_index, local_num_elems*sizeof(int));
    start_from = local_num_elems;
    if(output_style == VTK_ALL || output_style == VTK_RECV_LST) {
        for (i=0; i<recv_cnt[neighbour]; ++i) {
            local_global_index_big[start_from + i] = recv_lst[neighbour][i];
        }
        start_from += recv_cnt[neighbour];
    }
    if(output_style == VTK_ALL || output_style == VTK_SEND_LST) {
        for (i=0; i<send_cnt[neighbour]; ++i) {
            local_global_index_big[start_from + i] = send_lst[neighbour][i];
        }
    }

    // Fill scalars
    // Cells which belongs to process 0 will have value 1 to distinguish it from other cells
    // and other processors are incremented by one!
    for(i=0; i<local_num_elems; ++i) {
        scalars[i] = myrank+1;
    }
    start_from = local_num_elems;
    if(output_style == VTK_ALL || output_style == VTK_RECV_LST) {
        for (i=0; i<recv_cnt[neighbour]; ++i) {
            scalars[start_from + i] = nghb_to_rank[neighbour] + 1;
        }
        start_from += recv_cnt[neighbour];
    }
    if(output_style == VTK_ALL || output_style == VTK_SEND_LST) {
        for (i=0; i<send_cnt[neighbour]; ++i) {
            scalars[start_from + i] = nghb_to_rank[neighbour] + 0.5;
        }
    }

    // Write VTK
    char szFileName[80];
    //TODO:externalize this string
    const char *kOutputDirectoryName = "./out/";

    //find base file name
    char *data_file = strrchr(file_in,'/')+1;
    //strip data file base name
    data_file = strndup(data_file, strchr(data_file, '.')-data_file);

    if(output_style == VTK_SEND_LST) {
        sprintf(szFileName, "%s%s.Send.rank%i.neighbour%d.vtk", kOutputDirectoryName, data_file, myrank, neighbour);
    } else if(output_style == VTK_RECV_LST) {
        sprintf(szFileName, "%s%s.Recv.rank%i.neighbour%d.vtk", kOutputDirectoryName, data_file, myrank, neighbour);
    } else {
        sprintf(szFileName, "%s%s.RecvandSend.rank%d.neighbour%i.vtk", kOutputDirectoryName, data_file, myrank, neighbour);
    }


    test_distribution(file_in, szFileName, local_global_index_big, local_num_elems_big, scalars);
}
// end_of_student_code-----------------------------------------------------------------------------------
