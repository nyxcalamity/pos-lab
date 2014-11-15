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
echo "[INFO] Executing oneread strategies for 4 tasks"
mpiexec -n 4 ./gccg ./data/drall.geo.bin classic oneread
mpiexec -n 4 ./gccg ./data/drall.geo.bin nodal oneread
mpiexec -n 4 ./gccg ./data/drall.geo.bin dual oneread
echo "[INFO] Executing allread strategies for 4 tasks"
mpiexec -n 4 ./gccg ./data/drall.geo.bin classic allread
mpiexec -n 4 ./gccg ./data/drall.geo.bin nodal allread
mpiexec -n 4 ./gccg ./data/drall.geo.bin dual allread
echo "---------------------------------------------------"
echo "[INFO] Executing oneread strategies for 12 tasks"
mpiexec -n 12 ./gccg ./data/cojack.geo.bin classic oneread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin nodal oneread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin dual oneread
echo "[INFO] Executing allread strategies for 12 tasks"
mpiexec -n 12 ./gccg ./data/cojack.geo.bin classic allread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin nodal allread
mpiexec -n 12 ./gccg ./data/cojack.geo.bin dual allread