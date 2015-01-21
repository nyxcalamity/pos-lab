#/bin/sh
# Submits all statistics collection jobs to the LoadLeveler

#save directories
SCRIPT_DIR=$(dirname $0)
BASE_DIR=$SCRIPT_DIR/..
WDIR=$(pwd)

#store code location
CODE=$1
if [ -z "$CODE" ]; then
    CODE=$BASE_DIR/code
else 
    CODE=$BASE_DIR/$CODE
fi

#make the code
cd $CODE
make

#navigate back and clean up stats
cd $WDIR
rm -rf $BASE_DIR/stats/*

#generate job files
#NUM_PROC=(2 3 4 5 6 7 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64)
NUM_PROC=(68 72 76 80 84 88 96 112 128)
for i in ${NUM_PROC[@]}; do
    $SCRIPT_DIR/job-generator.sh $i $CODE $SCRIPT_DIR
done

#submit all jobs
for i in ${NUM_PROC[@]}; do 
    llsubmit $SCRIPT_DIR/a24-np$i.sh
done

#wait until jobs are finished
while [ ! -f $BASE_DIR/stats/job24-np128.out ]
do
  sleep 5
done

#pause for a bit and collect all the output
sleep 5
cat $BASE_DIR/stats/*.err >> $BASE_DIR/stats/all.err
cat $BASE_DIR/stats/*.out >> $BASE_DIR/stats/all.out