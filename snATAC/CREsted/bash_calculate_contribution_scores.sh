#!/bin/bash
#SBATCH --job-name=contrib_calcs
#SBATCH --output=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/analysis/log_files/modisco_ft_500_%j.out
#SBATCH --error=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/analysis/log_files/modisco_ft_500_%j.err
#SBATCH --partition=celltypes
#SBATCH --gpus=a100:1
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00


# === ENVIRONMENT SETUP ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate crested

# -------------- CONFIGURABLE PARAMS -----------------

## Human
adata_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/crested_adata/human_spinalcord_hmba_crested_chromsplit_norm.h5ad"
genome_file="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/fasta/genome.fa"
chrom_sizes="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/star/chrNameLength.txt"
model_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/analysis/spc_human/ft/checkpoints/06.keras"
output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/human/analysis/modisco/modisco_results_ft_500/"

## Macaque
#adata_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/macaque/crested_adata/macaque_spinalcord_hmba_crested_chromsplit_norm2.h5ad"
#genome_file="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/macaque/ncbi/mmul10/genome/fasta/genome.fa"
#chrom_sizes="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/macaque/ncbi/mmul10/genome/star/chrNameLength.txt"
#model_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/macaque/analysis/spc_macaque/ft3/checkpoints/05.keras"
#output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/macaque/analysis/modisco/modisco_results_ft_500_2_prediction/"

## Mouse
#adata_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/mouse/crested_adata/mouse_spinalcord_hmba_crested_chromsplit_norm.h5ad"
#genome_file="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/mouse/10x/mm10/genome/fasta/genome.fa"
#chrom_sizes="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/mouse/10x/mm10/genome/star/chrNameLength.txt"
#model_path="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/mouse/analysis/spc_mouse/ft2/checkpoints/24.keras"
#output_dir="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/spinal_cord/data/mouse/analysis/modisco/modisco_results_ft_500_prediction"

top_k=500
gradient_method='integrated_grad'
filtering_method='combined' # prediction, atac, or combined (default)

# -------------- LOG PARAMS TO STDOUT ----------------
echo "Starting run with the following parameters:"
echo "adata_path: $adata_path"
echo "genome_file: $genome_file"
echo "chrom_sizes: $chrom_sizes"
echo "model_path: $model_path"
echo "output_dir: $output_dir"
echo "top_k: $top_k"
echo "gradient_method : $gradient_method"
echo "filtering_method : $filtering_method"
echo "SLURM_JOB_ID: $SLURM_JOB_ID"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "--------------------------------------------------"

# -------------- RUN THE SCRIPT ----------------------
python calculate_contribution_scores.py \
  --adata_path "$adata_path" \
  --genome_file "$genome_file" \
  --chrom_sizes "$chrom_sizes" \
  --model_path "$model_path" \
  --output_dir "$output_dir" \
  --top_k "$top_k" \
  --gradient_method "$gradient_method" \
  --filtering_method "$filtering_method" \