#!/bin/bash
#@ job_name = pos_ws21
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
##@ total_tasks=2
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws21_tag

#@ initialdir = $(home)/pos-lab-den/A2.1/code
#@ output = $(home)/job.out
#@ error = $(home)/job.err

#@ notification=always
#@ notify_user=denys.sobchyshak@tum.de
#@ queue

#system setup
perf_off

##execution
echo "============================================================================================="
echo "Testing performance"
echo "============================================================================================="
echo "[INFO] Executing read strategies for 4 tasks and DRALL data"
mpiexec -n 4 ./gccg ./data/drall.geo.bin classic oneread
mpiexec -n 4 ./gccg ./data/drall.geo.bin nodal oneread
mpiexec -n 4 ./gccg ./data/drall.geo.bin dual oneread
mpiexec -n 4 ./gccg ./data/drall.geo.bin classic allread
mpiexec -n 4 ./gccg ./data/drall.geo.bin nodal allread
mpiexec -n 4 ./gccg ./data/drall.geo.bin dual allread
echo "---------------------------------------------------"
echo "[INFO] Executing read strategies for 12 tasks and COJACK data"
mpiexec -n 12 ./gccg ./data/cojack.geo.bin classic oneread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin nodal oneread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin dual oneread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin classic allread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin nodal allread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin dual allread
#echo "============================================================================================="
#echo "Testing partitioning"
#echo "============================================================================================="
#echo "[INFO] Executing strategies for 9 tasks and PENT data"
#mpiexec -n 9 ./gccg ./data/pent.geo.bin and classic allread
#mpiexec -n 9 ./gccg ./data/pent.geo.bin and nodal allread
#mpiexec -n 9 ./gccg ./data/pent.geo.bin and dual allread
#echo "---------------------------------------------------"
#echo "[INFO] Executing strategies for 9 tasks and PENT data"
#mpiexec -n 9 ./gccg ./data/cojack.geo.bin classic allread
#mpiexec -n 9 ./gccg ./data/cojack.geo.bin nodal allread
#mpiexec -n 9 ./gccg ./data/cojack.geo.bin dual allread
#echo "============================================================================================="