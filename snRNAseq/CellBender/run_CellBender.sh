#!/bin/bash
#SBATCH --job-name=cellbender_bg
#SBATCH --partition=celltypes 
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=3:30:00
#SBATCH --array=1-204%6
#SBATCH --output=logs/cellbender_%A_%a.out
#SBATCH --error=logs/cellbender_%A_%a.err

set -eo pipefail

pwd; hostname; date

source /home/yuanyuan.fu/.bashrc
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

conda activate cellbender


BASE=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Yuanyuan/BG_cellbender/aibs
RAW_DIR="${BASE}/01.rawdata"
OUT_DIR="${BASE}/02.removeBackground"
LIB_LIST="${BASE}/lib.txt"

lib=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIB_LIST" | tr -d '\r' | xargs)

echo "Running library: $lib"

mkdir -p "${OUT_DIR}/${lib}"

input=$(ls "${RAW_DIR}/${lib}_raw_feature_bc_matrix.h5" | head -n 1)

echo "Input: $input"
echo "Output dir: ${OUT_DIR}/${lib}"

cd "${OUT_DIR}/${lib}"

cellbender remove-background \
    --cuda \
    --input "$input" \
    --output "${lib}.out.h5" \
    > "${lib}.cellbender.log" 2>&1


