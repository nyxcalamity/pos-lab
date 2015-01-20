#!/bin/bash
#@ job_name = pos_ws24
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks = 12
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws24_tag

#@ initialdir = $(home)/pos-lab/A2.4
#@ output = $(home)/pos-lab/A2.4/stats/job24-np12.out
#@ error = $(home)/pos-lab/A2.4/stats/job24-np12.err

#@ notification=always
#@ notify_user=denys.sobchyshak@tum.de
#@ queue

#system setup
perf_off

#minimum optimization objective command
echo "[INFO] Statistics collection on 12 processes"
echo "----------------------------------------------------------------------------------------------"
echo "[INFO] Running base code"
mpiexec -n 12 ./code/gccg ./code/data/pent.geo.bin nodal oneread
mpiexec -n 12 ./code/gccg ./code/data/pent.geo.bin nodal oneread
mpiexec -n 12 ./code/gccg ./code/data/pent.geo.bin nodal oneread
