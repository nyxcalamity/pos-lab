#/bin/sh
# Submits all statistics collection jobs to the LoadLeveler

llsubmit $(dirname $0)/job-run-np2.sh
llsubmit $(dirname $0)/job-run-np3.sh
llsubmit $(dirname $0)/job-run-np4.sh
llsubmit $(dirname $0)/job-run-np5.sh
llsubmit $(dirname $0)/job-run-np6.sh
llsubmit $(dirname $0)/job-run-np7.sh
llsubmit $(dirname $0)/job-run-np8.sh