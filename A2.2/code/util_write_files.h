#ifndef XWRITE_H_
#define XWRITE_H_

/**
 * Export the simulation statistics to a text file
 *
 * @param in_file_name input file used for this simulation
 * @param out_file_name target output file
 * @param nintci index of the first element to include in the output file
 * @param nintcf index of the last element to include in the output file
 * @param var array with the computed results
 * @param total_iters total number of iterations required to compute the solution
 * @param residual_ratio the final residual ratio
 * @return
 */
int store_simulation_stats(char *in_file_name, char *out_file_name, int nintci, int nintcf,
                           double *var, int total_iters, double residual_ratio);

/**
 * Write a VTK header with the desired geometry data to a file
 *
 * @param experiment_name title of the VTK file (in out case simply the used input file)
 * @param out_file_name target output file
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param points_count total number of points
 * @param points coordinates of the points that define the volume cells
 * @param elems neighboring information for the cells
 */
void vtk_write_unstr_grid_header(const char *experiment_name, const char *out_file_name,
                                 int start_index, int end_index, int points_count, int **points,
                                 int *elems);

/**
 * Append the values of a double variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
void vtk_append_double(const char *out_file_name, const char *var_name, int start_index,
                       int end_index, double *values);

/**
 * Append the values of an integer variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
void vtk_append_integer(const char *out_file_name, const char *var_name, int start_index,
                        int end_index, int *values);

// start_of_student_code---------------------------------------------------------------------------------
/**
 * Create VTK output files using data which is only on MYRANK processor
 *
 * @param experiment_name title of the VTK file (in out case simply the used input file)
 * @param out_file_name target output file
 * @param nintci_l index of the first element to include in the output file
 * @param nintcf_l index of the last element to include in the output file
 * @param points_count total number of points
 * @param points coordinates of the points that define the volume cells
 * @param elems neighboring information for the cells
 */
void vtk_for_process(const char *file_in, const char *file_vtk_out,
        int nintci, int nintcf, int points_count, int **points,  int *elems,
		int *local_global_index, int local_num_elems, double *scalars);


int vtk_check(char *file_in, int myrank,
        int nintci, int nintcf, double *su, double *cgup,
		int points_count, int **points,  int *elems,
		int *local_global_index, int local_num_elems);


// FIXME: comments
void vtk_check_lists(char *file_in, int myrank,
        int *local_global_index, int local_num_elems,
        int nghb_cnt, int* nghb_to_rank, int* send_cnt, int** send_lst,
        int *recv_cnt, int** recv_lst, int output_style);

void vtk_check_neighbour(char *file_in, int myrank,
        int *local_global_index, int local_num_elems,
        int nghb_cnt, int* nghb_to_rank, int* send_cnt, int** send_lst,
        int *recv_cnt, int** recv_lst, int output_style, int neighbour);
// end_of_student_code-----------------------------------------------------------------------------------
#endif /* XWRITE_H_ */

