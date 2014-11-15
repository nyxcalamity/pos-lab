####################################################################
#             Makefile for FIRE Solver Benchmarks                  #
####################################################################
CC = mpicc
CFLAGS = -g -O0 -Wall
LIBS = -lm -lmetis -lpapi

CFLAGS += $(METIS_INC)
LIBS += $(METIS_LIB)

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c \
       initialization_algorithms.c util_processors.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: clean init gccg

init:
	mkdir -p ./out

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c test_functions.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg $(LIBPOS) out/* pstats.dat
