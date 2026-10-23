#!/bin/bash
# Run guidance/01_reduce_fiber_field.py on one field per array task.
# guidance/fields_todo.txt (not in the repo) has one "<category> <field>" per line.
# Usage, from the repo root: sbatch --array=1-$(wc -l < guidance/fields_todo.txt)%6 guidance/run_reduce.sbatch
# Keep it at %6: with 12 concurrent tasks they starve each other on I/O.
#SBATCH --job-name=ecmreduce
#SBATCH --partition=general
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
# one 950 MB pixel CSV is held in memory at a time
#SBATCH --time=10:00:00
#SBATCH --output=logs/slurm_reduce_%A_%a.out

ROOT=.
PY=python
LIST=$ROOT/guidance/fields_todo.txt

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)
CATEGORY=$(echo $LINE | awk '{print $1}')
FIELD=$(echo $LINE | awk '{print $2}')

echo "task $SLURM_ARRAY_TASK_ID  node $(hostname)  $CATEGORY $FIELD  start $(date)"
$PY $ROOT/guidance/01_reduce_fiber_field.py $CATEGORY $FIELD
echo "exit=$?  end $(date)"
