#ifndef _POSL_DEFINITIONS_H_
#define _POSL_DEFINITIONS_H_

//quantitative definitions of cmd params
#define POSL_DATA_TJUNC             1
#define POSL_DATA_DRALL             2
#define POSL_DATA_PENT              3
#define POSL_DATA_COJACK            4

#define POSL_PARTITIONING_CLASSIC   1
#define POSL_PARTITIONING_DUAL      2
#define POSL_PARTITIONING_NODAL     3

#define POSL_INIT_ONE_READ          1
#define POSL_INIT_ALL_READ          2

//mpi communication tags
#define POSL_MPI_TAG_NINTCF         0
#define POSL_MPI_TAG_CELLS          200
#define POSL_MPI_TAG_NEXTCF         1
#define POSL_MPI_TAG_POINTS_COUNT   2
#define POSL_MPI_TAG_LCC            3
#define POSL_MPI_TAG_VAR            4
#define POSL_MPI_TAG_L2G            5
#define POSL_MPI_TAG_POINTS         99
#define POSL_MPI_TAG_ELEMENTS       100
#define POSL_MPI_TAG_BS             101
#define POSL_MPI_TAG_BE             102
#define POSL_MPI_TAG_BN             103
#define POSL_MPI_TAG_BW             104
#define POSL_MPI_TAG_BL             105
#define POSL_MPI_TAG_BH             106
#define POSL_MPI_TAG_BP             107
#define POSL_MPI_TAG_SU             108

//statistics output
#define POSL_COLLECT_STATS          1

//error definitions
#define POSL_ERROR                  -1
#define POSL_OK                     0

//debug definitions
#define DEBUG_ENABLED               0
#define OUTPUT_LCC_G                0
#define OUTPUT_LCC                  0
#define OUTPUT_NINTCF_NINTCE        0
#define DEBUG_OUTPUT_L2G_G          0
#define DEBUG_OUTPUT_PARTITIONING   0
        
// VTK output definitions
#define POSL_OUTPUT_DIRECTORY       "./out/"
#define POSL_OUTPUT_DIR             "out"
#define OUTPUT_VTK                  0
#define VTK_ALL                     1
#define VTK_SEND_LST                2
#define VTK_RECV_LST                3

#define VTK_NEIGHBOUR               1

#endif
