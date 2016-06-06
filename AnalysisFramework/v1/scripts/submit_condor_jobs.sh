#!/bin/bash
# [1]==first run number
# [2]==last run number

#set coredump file sizes to zero
ulimit -c 0

for (( i=$1; i<=$2; i++ ))
do
  export filename=/net/mitbbr00/data03/ddujmic/data/dmtpc_run00$i.root;
  if [ -f $filename ];
  then
     echo "Submitting condor job: $i ..." $filename;
     ./make_condor_job.py $i;
     condor_submit -d 'requirements = Machine != "qweak.lns.mit.edu"' job_scripts/condor_spec_$i;
  fi
done



