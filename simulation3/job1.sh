#!/bin/bash
#SBATCH -J p06s3j1
#SBATCH  --nodes=1
#SBATCH  --ntasks=7
#SBATCH --cpus-per-task=1
#SBATCH --time=720:00:00

#SBATCH --mail-user=zhubencong@gmail.com
#SBATCH -o simulation3/jobout/job1_out.log
#SBATCH -e simulation3/jobout/job1_err.log


source /opt/share/etc/miniconda3-py39.sh
conda activate R4.2
export OMP_NUM_THREADS=1


cd /home/project07/Bencong/scratch/project06


for p in 2000 5000 10000 15000 20000 25000 30000
do
    srun --ntasks=1 --exclusive R --vanilla --slave --args $p  < simulation3/sim3DataGenerator2.R & 
done
wait