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
#@ output = $(home)/job-papi.out
#@ error = $(home)/job-papi.err

#@ notification=always
#@ notify_user=denys.sobchyshak@tum.de
#@ queue

#system setup
perf_off

##execution
echo "-------------------------------------------------------------------------"
echo "Testing performance with -O1 flag."
echo "-------------------------------------------------------------------------"
echo "[INFO]Running tests with cojack.dat. Trial 1."
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/cojack.dat -cache
echo "[INFO]Running tests with cojack.dat. Trial 2."
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/cojack.dat -cache
echo "[INFO]Running tests with cojack.dat. Trial 3."
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/cojack.dat -mflops
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/cojack.dat -cache
echo "-------------------------------------------------------------------------"
echo "[INFO]Running tests with tjunc.dat. Trial 1."
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/tjunc.dat -mflops
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/tjunc.dat -cache
echo "[INFO]Running tests with tjunc.dat. Trial 2."
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/tjunc.dat -mflops
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/tjunc.dat -cache
echo "[INFO]Running tests with tjunc.dat. Trial 3."
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/tjunc.dat -mflops
pos-lab/A1/code-papi-o1/gccg pos-lab/A1/code-papi-o1/data/tjunc.dat -cache

echo "-------------------------------------------------------------------------"
echo "Testing performance with -O3 flag."
echo "-------------------------------------------------------------------------"
echo "[INFO]Running tests with cojack.dat. Trial 1."
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/cojack.dat -cache
echo "[INFO]Running tests with cojack.dat. Trial 2."
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/cojack.dat -cache
echo "[INFO]Running tests with cojack.dat. Trial 3."
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/cojack.dat -mflops
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/cojack.dat -cache
echo "-------------------------------------------------------------------------"
echo "[INFO]Running tests with tjunc.dat. Trial 1."
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/tjunc.dat -mflops
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/tjunc.dat -cache
echo "[INFO]Running tests with tjunc.dat. Trial 2."
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/tjunc.dat -mflops
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/tjunc.dat -cache
echo "[INFO]Running tests with tjunc.dat. Trial 3."
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/tjunc.dat -mflops
pos-lab/A1/code-papi-o3/gccg pos-lab/A1/code-papi-o3/data/tjunc.dat -cache