#/bin/sh
# Submits all statistics collection jobs to the LoadLeveler

# Clean all stats dirs
rm -rf $(dirname $0)/../stats/*

#submit all jobs
llsubmit $(dirname $0)/job-run-np12.sh
llsubmit $(dirname $0)/job-run-np16.sh
llsubmit $(dirname $0)/job-run-np20.sh
llsubmit $(dirname $0)/job-run-np24.sh
llsubmit $(dirname $0)/job-run-np28.sh
llsubmit $(dirname $0)/job-run-np32.sh
llsubmit $(dirname $0)/job-run-np36.sh

#wait until jobs are finished
while [ ! -f $(dirname $0)/../stats/job24-np36.out ]
do
  sleep 2
done

#pause for a bit and collect all the output
sleep 4
cat $(dirname $0)/../stats/*.err >> $(dirname $0)/../stats/all.err
cat $(dirname $0)/../stats/*.out >> $(dirname $0)/../stats/all.out
