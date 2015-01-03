#!/bin/bash
#@ job_name = pos_ws24
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks = 9
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws24_tag

#@ initialdir = $(home)/pos-lab/A2.4
#@ output = $(home)/job24.out
#@ error = $(home)/job24.err

#@ notification=always
#@ notify_user=denys.sobchyshak@tum.de
#@ queue

#system setup
perf_off

#minimum optimization objective command
mpiexec -n 9 ./optimization1/optimized/gccg ./code/data/pent.geo.bin nodal oneread

