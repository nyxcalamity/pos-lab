# Makefile for FIRE Solver Benchmark
#===================================================================================================
# Declaration of variables
#===================================================================================================
CC = mpicc
CFLAGS = -Wall -O3 #-g
LIBS = -lm

# Dynamically add
CFLAGS += $(METIS_INC) $(PAPI_INC)
LIBS += $(METIS_LIB) $(PAPI_LIB)

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c \
       initialization_algorithms.c util_processors.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

#===================================================================================================
# Cleanup and initialization tasks
#===================================================================================================
all: clean init gccg

init:
	mkdir -p ./out
	
clean: clean-data
	rm -rf *.o gccg $(LIBPOS)

clean-data:
	rm -rf out/* *.dat data/*.out

#===================================================================================================
# Compilation tasks
#===================================================================================================
%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c test_functions.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+