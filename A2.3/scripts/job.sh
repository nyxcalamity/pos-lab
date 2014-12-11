#!/bin/bash
#@ job_name = pos_ws22
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks = 9
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws22_tag

#@ initialdir = $(home)/pos-lab/A2.3/code
#@ output = $(home)/job.out
#@ error = $(home)/job.err

#@ notification=always
#@ notify_user=denys.korzh@tum.de
#@ queue

#system setup
perf_off

##execution
echo "============================================================================================="
echo "Statistics output"
echo "============================================================================================="
echo "============================================================================================="
echo "drall output"
echo "============================================================================================="
mpiexec -n 9 ./gccg ./data/drall.geo.bin classic allread
mpiexec -n 9 ./gccg ./data/drall.geo.bin classic oneread
mpiexec -n 9 ./gccg ./data/drall.geo.bin dual allread
mpiexec -n 9 ./gccg ./data/drall.geo.bin dual oneread
mpiexec -n 9 ./gccg ./data/drall.geo.bin nodal allread
mpiexec -n 9 ./gccg ./data/drall.geo.bin nodal oneread
echo "============================================================================================="
echo "pent output"
echo "============================================================================================="
mpiexec -n 9 ./gccg ./data/pent.geo.bin classic allread
mpiexec -n 9 ./gccg ./data/pent.geo.bin classic oneread
mpiexec -n 9 ./gccg ./data/pent.geo.bin dual allread
mpiexec -n 9 ./gccg ./data/pent.geo.bin dual oneread
mpiexec -n 9 ./gccg ./data/pent.geo.bin nodal allread
mpiexec -n 9 ./gccg ./data/pent.geo.bin nodal oneread
echo "============================================================================================="
echo "cojack output"
echo "============================================================================================="
mpiexec -n 9 ./gccg ./data/cojack.geo.bin classic allread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin classic oneread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin dual allread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin dual oneread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin nodal allread
mpiexec -n 9 ./gccg ./data/cojack.geo.bin nodal oneread
echo "============================================================================================="
echo "Finish output"
echo "============================================================================================="
