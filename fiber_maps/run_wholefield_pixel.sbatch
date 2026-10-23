#!/bin/bash
#SBATCH --job-name=wf_pixel
#SBATCH --partition=general
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
# 32G is too tight: a 48-frame 4923x4927 field peaked at 32.0 GiB (memory grows over the run)
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm_wfpixel_%j.out

# Run alignment_analysis_pixel.py on all frames of one whole field at native resolution.
# usage (from the repo root): sbatch fiber_maps/run_wholefield_pixel.sbatch P201015 02
# Raw tifs are symlinked into results/fib_wholefield/<category>/<sample>_field_<NN>/,
# and all UMA output is written there.

SAMPLE=$1
FIELD=$2

if [ -z "$SAMPLE" ] || [ -z "$FIELD" ]; then
    echo "usage: sbatch run_wholefield_pixel.sbatch <SAMPLE> <FIELD>   e.g. P201015 02"
    exit 1
fi

ROOT=.
PY=python
SCRIPT=$ROOT/fiber_maps/alignment_analysis_pixel.py
NAME=${SAMPLE}_field_${FIELD}

# locate the raw fibronectin folder (it sits under one of the 4 category folders)
SRC=$(ls -d $ROOT/data/data_raw/data_raw_fib/*/$NAME 2>/dev/null)
if [ $(echo "$SRC" | wc -l) -ne 1 ] || [ ! -d "$SRC" ]; then
    echo "ERROR: expected exactly one $NAME under data_raw_fib/*/, got: $SRC"
    exit 1
fi
CATEGORY=$(basename $(dirname $SRC))

WORK=$ROOT/results/fib_wholefield/$CATEGORY/$NAME
mkdir -p $WORK
ln -sfn $SRC/*.tif $WORK/

cat > $WORK/input_paths.json <<EOF
{
  "folder_paths": ["$WORK"]
}
EOF

echo "node:     $(hostname)"
echo "sample:   $NAME  (category: $CATEGORY)"
echo "source:   $SRC"
echo "frames:   $(ls $SRC/*.tif | wc -l)"
echo "work dir: $WORK"
echo "start:    $(date)"
echo "-----------------------------------------------------------"

# '1' = fibronectin channel index (these tifs are single channel), 'y' = start
printf '1\ny\n' | $PY $SCRIPT -i $WORK/input_paths.json

echo "-----------------------------------------------------------"
echo "exit=$?  end: $(date)"
