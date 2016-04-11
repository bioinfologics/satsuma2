# Script for starting Satsuma jobs on different job submission environments
# Comment out the lines not required
# Usage: satsuma_run.sh <current_path> <kmatch_cmd> <ncpus> <mem> <job_id> <run_synchronously>
# mem should be in Gb, ie. 100Gb = 100

# no submission system, run process locally either synchronously or asynchronously
if [ "$6" -eq 1 ]; then
  eval "$2"
else
  eval "$2" &
fi

# qsub (PBS systems)
#echo "cd $1; $2" | qsub -V -qProd -l ncpus=$3,mem=$4G -N $5

# bsub (LSF systems)
#mem=`expr $4 + 1000`
#bsub -o ${5}.log -J $5 -n $3 -q Prod256 -R "rusage[mem=$mem]" "$2"

# SLURM systems
#echo "#!/bin/sh" > slurm_tmp.sh
#echo srun $2 >> slurm_tmp.sh
#sbatch -p tgac-long -c $3 -J $5 -o ${5}.log --mem ${4}G slurm_tmp.sh
