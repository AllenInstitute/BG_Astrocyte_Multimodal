#!/bin/bash
#SBATCH --job-name=run_modisco
#SBATCH --output=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/macaque/analysis/modisco/run_modisco_intgrad_%j.out
#SBATCH --error=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/macaque/analysis/modisco/run_modisco_intgrad_%j.err
#SBATCH --partition=celltypes
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=36:00:00

# === ENVIRONMENT SETUP ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate crested

output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/analysis/modisco/modisco_results_ft_500/"
contrib_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/analysis/modisco/modisco_results_ft_500/"
window=1000
fdr=0.2
max_seqlets=20000
report=false

# -------------- LOG PARAMS TO STDOUT ----------------
echo "Starting run with the following parameters:"
echo "output_dir: $output_dir"
echo "contrib_dir: $contrib_dir"
echo "window: $window"
echo "fdr: $fdr"
echo "max_seqlets: $max_seqlets"
echo "report: $report"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "--------------------------------------------------"

# -------------------- RUN SCRIPT --------------------
python /home/niklas.kempynck/nkemp/spc_crested/run_modisco.py \
    --window $window \
    --output_dir $output_dir \
    --contrib_dir $contrib_dir \
    --fdr $fdr \
    --max_seqlets $max_seqlets \
    --report
