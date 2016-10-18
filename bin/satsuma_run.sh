#!/bin/sh

# Script for starting Satsuma jobs on different job submission environments
# One section only should be active, ie. not commented out

# Usage: satsuma_run.sh <current_path> <kmatch_cmd> <ncpus> <mem> <job_id> <run_synchronously>
# mem should be in Gb, ie. 100Gb = 100

# no submission system, processes are run locally either synchronously or asynchronously
if [ "$6" -eq 1 ]; then
  eval "$2"
else
  eval "$2" &
fi

##############################################################################################################
## For the sections below you will need to change the queue name (QueueName) to one existing on your system ##
##############################################################################################################

# qsub (PBS systems)
#echo "cd $1; $2" | qsub -V -qQueueName -l ncpus=$3,mem=$4G -N $5

# bsub (LSF systems)
#mem=`expr $4 + 1000`
#bsub -o ${5}.log -J $5 -n $3 -q QueueName -R "rusage[mem=$mem]" "$2"

# SLURM systems
#echo "#!/bin/sh" > slurm_tmp.sh
#echo srun $2 >> slurm_tmp.sh
#sbatch -p QueueName -c $3 -J $5 -o ${5}.log --mem ${4}G slurm_tmp.sh
