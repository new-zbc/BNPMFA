#!/bin/bash
#SBATCH -J p06s3j3
#SBATCH  --nodes=1
#SBATCH  --ntasks=7
#SBATCH --cpus-per-task=1
#SBATCH --time=720:00:00

#SBATCH --mail-user=zhubencong@gmail.com
#SBATCH -o simulation3/jobout/job3_out.log
#SBATCH -e simulation3/jobout/job3_err.log


source /opt/share/etc/miniconda3-py39.sh
conda activate R4.2
export OMP_NUM_THREADS=1


cd /home/project07/Bencong/scratch/project06


for DataFolder in "p_2000" "p_5000" "p_10000" "p_15000" "p_20000" "p_25000" "p_30000"
do
  # srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder  < simulation3/script/BayesSpace.R & 
  # srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder  < simulation3/script/BANKSY.R & 
  # srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder  < simulation3/script/DRSC.R & 
  # srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder  < simulation3/script/SCMEB.R & 
  # srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder  < simulation3/script/construct_h5ad.R & 
  srun --ntasks=1 --exclusive R --vanilla --slave --args $DataFolder  < simulation3/script/Louvain.R & 
  
done
wait