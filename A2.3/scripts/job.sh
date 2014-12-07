#!/bin/bash
#@ job_name = pos_ws22
#@ job_type = MPICH
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks = 9
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws22_tag

#@ initialdir = $(home)/pos-lab/A2.2/code
#@ output = $(home)/job.out
#@ error = $(home)/job.err

#@ notification=always
#@ notify_user=denys.sobchyshak@tum.de
#@ queue

#system setup
perf_off

##execution
echo "============================================================================================="
echo "Statistics output"
echo "============================================================================================="
mpiexec -n 9 ./gccg ./data/drall.geo.bin classic allread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin classic allread

mpiexec -n 9 ./gccg ./data/drall.geo.bin dual oneread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin dual oneread