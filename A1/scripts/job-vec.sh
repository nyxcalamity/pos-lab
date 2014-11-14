#!/bin/bash
#@ job_name = pos_ws1_papi
#@ job_type = parallel
#@ class = test
#@ island_count=1
#@ node = 1
#@ total_tasks=16
#@ wall_clock_limit = 0:30:00
#@ energy_policy_tag = pos_ws1_papi_tag
##@ minimize_time_to_solution = yes

#@ initialdir = $(home)
#@ output = $(home)/job-vec.out
#@ error = $(home)/job-vec.err

#@ notification=always
#@ notify_user=denys.sobchyshak@tum.de
#@ queue

#system setup
perf_off

##execution
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O1 -no-vec"
echo "-------------------------------------------------------------------------"
pos-lab/A1/code-papi-o1-novec/gccg pos-lab/A1/code-papi-o1-novec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1-novec/gccg pos-lab/A1/code-papi-o1-novec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1-novec/gccg pos-lab/A1/code-papi-o1-novec/data/cojack.dat -mflops
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O1 -vec"
echo "-------------------------------------------------------------------------"
pos-lab/A1/code-papi-o1-vec/gccg pos-lab/A1/code-papi-o1-vec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1-vec/gccg pos-lab/A1/code-papi-o1-vec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1-vec/gccg pos-lab/A1/code-papi-o1-vec/data/cojack.dat -mflops
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O1 -vec -xhost"
echo "-------------------------------------------------------------------------"
pos-lab/A1/code-papi-o1-vec-xhost/gccg pos-lab/A1/code-papi-o1-vec-xhost/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1-vec-xhost/gccg pos-lab/A1/code-papi-o1-vec-xhost/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1-vec-xhost/gccg pos-lab/A1/code-papi-o1-vec-xhost/data/cojack.dat -mflops
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O3 -no-vec"
echo "-------------------------------------------------------------------------"
pos-lab/A1/code-papi-o3-novec/gccg pos-lab/A1/code-papi-o3-novec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3-novec/gccg pos-lab/A1/code-papi-o3-novec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3-novec/gccg pos-lab/A1/code-papi-o3-novec/data/cojack.dat -mflops
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O3 -vec"
echo "-------------------------------------------------------------------------"
pos-lab/A1/code-papi-o3-vec/gccg pos-lab/A1/code-papi-o3-vec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3-vec/gccg pos-lab/A1/code-papi-o3-vec/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3-vec/gccg pos-lab/A1/code-papi-o3-vec/data/cojack.dat -mflops
echo "-------------------------------------------------------------------------"
echo "Vectorization tests for optimization flags: -O3 -vec -xhost"
echo "-------------------------------------------------------------------------"
pos-lab/A1/code-papi-o3-vec-xhost/gccg pos-lab/A1/code-papi-o3-vec-xhost/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3-vec-xhost/gccg pos-lab/A1/code-papi-o3-vec-xhost/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3-vec-xhost/gccg pos-lab/A1/code-papi-o3-vec-xhost/data/cojack.dat -mflops