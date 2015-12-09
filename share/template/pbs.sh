#!/bin/tcsh
#PBS -N TEMPLATE_JOB_NAME
#PBS -o TEMPLATE_NAME.out
#PBS -j oe
#PBS -q TEMPLATE_PBS_QUEUE
source /ihepbatch/bes/nikolaev/bin/boss664
setenv exe 'boss.exe' 

set JOB_NAME=TEMPLATE_NAME
cd  TEMPLATE_DIR
set SIM_OPTS=$JOB_NAME"-sim.cfg"
set REC_OPTS=$JOB_NAME"-rec.cfg"
set SEL_OPTS=$JOB_NAME"-sel.cfg"
set LOG=$JOB_NAME.out
set ERR_LOG=$JOB_NAME.err

if ( -f $LOG ) then
  mv $LOG $LOG"_bak" 
endif
if ( -f $ERR_LOG ) then
   mv $ERR_LOG $ERR_LOG"_bak" 
endif

#((time ) >> $LOG) >& $ERR_LOG
#((time ) >> $LOG) >& $ERR_LOG
#((time ) >> $LOG) >& $ERR_LOG
$exe $SIM_OPTS
$exe $REC_OPTS
$exe $SEL_OPTS
