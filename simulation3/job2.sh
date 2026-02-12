#!/bin/bash
#SBATCH -J p06s3j2
#SBATCH  --nodes=1
#SBATCH  --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --time=720:00:00

#SBATCH --mail-user=zhubencong@gmail.com
#SBATCH -o simulation3/jobout/job2_out.log
#SBATCH -e simulation3/jobout/job2_err.log


source /opt/share/etc/miniconda3-py39.sh
conda activate R4.2
export OMP_NUM_THREADS=1


cd /home/project07/Bencong/scratch/project06


for DataFolder in "p_2000" "p_5000" "p_10000" "p_15000" "p_20000" "p_25000" "p_30000"
do
  for method in "PY" "DP" "MFM"
  do
    for f in 0 0.5 1 1.5 2 2.5 3 3.5
    do
      srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder $method $f  < simulation3/DRMFMsim.R & 
    done
  done
  wait
  echo "$DataFolder done"
done
wait