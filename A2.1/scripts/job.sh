#!/bin/bash
#@ job_name = pos_ws21
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks=2
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
mpiexec ./gccg ./data/pent.geo.bin nodal oneread