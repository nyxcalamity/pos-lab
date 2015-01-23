#!/bin/bash
#@ job_name = pos_ws21_seq
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks = 2
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws21_seq

#@ initialdir = $(home)/pos-lab/A2.1-sequential/code
#@ output = $(home)/job.out
#@ error = $(home)/job.err

#@ notification=always
#@ notify_user=denys.korzh@tum.de
#@ queue

#system setup
perf_off

##execution
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O3"
echo "-------------------------------------------------------------------------"
echo "PENT"
mpiexec -n 1 ./gccg ./data/pent.geo.bin classic oneread
mpiexec -n 1 ./gccg ./data/pent.geo.bin classic oneread
mpiexec -n 1 ./gccg ./data/pent.geo.bin classic oneread
echo "PENT"
mpiexec -n 1 ./gccg ./data/drall.geo.bin classic allread
mpiexec -n 1 ./gccg ./data/drall.geo.bin classic allread
mpiexec -n 1 ./gccg ./data/drall.geo.bin classic allread
echo "PENT"
mpiexec -n 1 ./gccg ./data/cojack.geo.bin classic allread
mpiexec -n 1 ./gccg ./data/cojack.geo.bin classic allread
mpiexec -n 1 ./gccg ./data/cojack.geo.bin classic allread