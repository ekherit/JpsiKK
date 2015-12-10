#!/bin/tcsh
#PBS -N TEMPLATE_JOB_NAME
#PBS -o TEMPLATE_NAME.out
#PBS -j oe
#PBS -q TEMPLATE_PBS_QUEUE
source /ihepbatch/bes/nikolaev/bin/boss664

set JOB_NAME=TEMPLATE_NAME

cd  TEMPLATE_DIR
