#/bin/sh
# Submits all statistics collection jobs to the LoadLeveler

# Clean all stats dirs
rm -rf $(dirname $0)/../stats/*

#submit all jobs
llsubmit $(dirname $0)/job-run-np2.sh
llsubmit $(dirname $0)/job-run-np3.sh
llsubmit $(dirname $0)/job-run-np4.sh
llsubmit $(dirname $0)/job-run-np5.sh
llsubmit $(dirname $0)/job-run-np6.sh
llsubmit $(dirname $0)/job-run-np7.sh
llsubmit $(dirname $0)/job-run-np8.sh

#wait until jobs are finished
while [ ! -f $(dirname $0)/../stats/job24-np8.out ]
do
  sleep 2
done

#pause for a bit and collect all the output
sleep 4
cat $(dirname $0)/../stats/*.err >> $(dirname $0)/../stats/all.err
cat $(dirname $0)/../stats/*.out >> $(dirname $0)/../stats/all.out