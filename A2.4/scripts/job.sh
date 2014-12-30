#!/bin/bash
#@ job_name = pos_ws23
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks = 9
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws23_tag

#@ initialdir = $(home)/pos-lab/A2.3/code
#@ output = $(home)/job.out
#@ error = $(home)/job.err

#@ notification=always
#@ notify_user=denys.korzh@tum.de
#@ queue

#system setup
perf_off

##execution
mpiexec -n 9 ./gccg ./data/drall.geo.bin classic oneread

