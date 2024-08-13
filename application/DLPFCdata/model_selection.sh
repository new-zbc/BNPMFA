#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -N p6-ms
#PBS -q batch
#PBS -o DLPFCdata/jobout/job1_out.log
#PBS -e DLPFCdata/jobout/job1_err.log

source /home/project07/anaconda3/etc/profile.d/conda.sh

conda activate R4.2

cd /home/project07/Bencong/project06

#R --vanilla --slave --args 151507 1 15 < DLPFCdata/model_selection.R &
#R --vanilla --slave --args 151508 3 7 < DLPFCdata/model_selection.R &

R --vanilla --slave --args 151509 1 15 < DLPFCdata/model_selection.R &
#R --vanilla --slave --args 151510 1 15 < DLPFCdata/model_selection.R &

#R --vanilla --slave --args 151669 2 15 < DLPFCdata/model_selection.R &
#R --vanilla --slave --args 151670 3 15 < DLPFCdata/model_selection.R &

#R --vanilla --slave --args 151671 3 15 < DLPFCdata/model_selection.R &
#R --vanilla --slave --args 151672 1 5 < DLPFCdata/model_selection.R &
wait
