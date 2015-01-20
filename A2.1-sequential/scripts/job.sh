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
echo "Vectorization tests for optimization flags: -O1 -no-vec"
echo "-------------------------------------------------------------------------"
mpiexec -n 1 ./gccg ./data/pent.geo.bin classic oneread


