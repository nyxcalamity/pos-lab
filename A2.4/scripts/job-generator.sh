#!/bin/bash
# Generates a job file with predefined number of tasks and nodes used

NP=$1
NN=$((($NP+15)/16)) 
EXEC=$2/gccg
FILE_NAME=$3/a24-np$NP.sh
rm -f $FILE_NAME

echo "#!/bin/bash" >> $FILE_NAME
echo "#@ job_name = pos_ws24" >> $FILE_NAME
echo "#@ job_type = parallel" >> $FILE_NAME
echo "#@ class = test" >> $FILE_NAME
echo "#@ island_count=1" >> $FILE_NAME
echo "#@ node = $NN" >> $FILE_NAME
echo "#@ total_tasks = $NP" >> $FILE_NAME
echo "#@ wall_clock_limit = 0:30:00" >> $FILE_NAME
echo "#@ energy_policy_tag = pos_ws24_tag" >> $FILE_NAME
echo "" >> $FILE_NAME
echo "#@ initialdir = \$(home)/pos-lab/A2.4" >> $FILE_NAME
echo "#@ output = \$(home)/pos-lab/A2.4/stats/job24-np$NP.out" >> $FILE_NAME
echo "#@ error = \$(home)/pos-lab/A2.4/stats/job24-np$NP.err" >> $FILE_NAME
echo "" >> $FILE_NAME
echo "#@ notification=always" >> $FILE_NAME
echo "#@ notify_user=denys.sobchyshak@tum.de" >> $FILE_NAME
echo "#@ queue" >> $FILE_NAME
echo "" >> $FILE_NAME
echo "#system setup" >> $FILE_NAME
echo "perf_off" >> $FILE_NAME
echo "" >> $FILE_NAME
echo "#minimum optimization objective command" >> $FILE_NAME
echo "echo \"[INFO] Statistics collection on $NP processes\"" >> $FILE_NAME
echo "echo \"---------------------------------------------------------------------------------------\"" >> $FILE_NAME
echo "echo \"Testing PENT data\"" >> $FILE_NAME
echo "echo \"---------------------------------------------------------------------------------------\"" >> $FILE_NAME
echo "echo \"[INFO] NODAL\"" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin nodal oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin nodal oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin nodal oneread" >> $FILE_NAME
echo "echo \"[INFO] DUAL\"" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin dual oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin dual oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin dual oneread" >> $FILE_NAME
echo "echo \"[INFO] CLASSIC\"" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin classic oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin classic oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/pent.geo.bin classic oneread" >> $FILE_NAME
echo "echo \"---------------------------------------------------------------------------------------\"" >> $FILE_NAME
echo "echo \"Testing COJACK data\"" >> $FILE_NAME
echo "echo \"---------------------------------------------------------------------------------------\"" >> $FILE_NAME
echo "echo \"[INFO] NODAL\"" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin nodal oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin nodal oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin nodal oneread" >> $FILE_NAME
echo "echo \"[INFO] DUAL\"" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin dual oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin dual oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin dual oneread" >> $FILE_NAME
echo "echo \"[INFO] CLASSIC\"" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin classic oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin classic oneread" >> $FILE_NAME
echo "mpiexec -n $NP $EXEC ./code/data/cojack.geo.bin classic oneread" >> $FILE_NAME