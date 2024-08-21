#!/bin/bash
#SBATCH -J p06-star
#SBATCH  --nodes=1
#SBATCH  --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --time=720:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=1155149265@link.cuhk.edu.hk
#SBATCH -o STARmap/jobout/ml_out.log
#SBATCH -e STARmap/jobout/ml_err.log

cd /home/project07/Bencong/project07
source /opt/share/etc/miniconda3-py39.sh
conda activate R4.2

cd /home/project07/Bencong/scratch/project06

#srun --ntasks=1 --exclusive R --vanilla --slave --args "BZ5" < STARmap/model_selection.R &
#srun --ntasks=1 --exclusive R --vanilla --slave --args "BZ9" < STARmap/model_selection.R &
srun --ntasks=1 --exclusive R --vanilla --slave --args "BZ14" < STARmap/model_selection.R &
wait

echo "out"
