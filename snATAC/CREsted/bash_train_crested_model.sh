#!/bin/bash
#SBATCH --job-name=train_crested
#SBATCH --output=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/log_files/cluster_log_train_crested_full_universe_subset_filtered_dars_%j.out
#SBATCH --error=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/log_files/cluster_log_train_crested_full_universe_subset_filtered_dars_%j.err
#SBATCH --partition=celltypes
#SBATCH --gpus=a100:1
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=48:00:00

# === CONFIGURATION ===
CONDA_ENV="crested_testing"

#Human
GENOME_FASTA="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/fasta/genome.fa"
CHR_SIZES="/allen/programs/celltypes/workgroups/rnaseqanalysis/references/human/10x/grch38.p2/genome/star/chrNameLength.txt"
ADATA_PATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/astro_cluster_subset_filtered.h5ad"
ADATA_PATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/dorsal_ventral_dars.h5ad"

PROJECT_NAME="astro_human"
cd /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/

RUN_NAME="cluster_full_universe_base_subset_filtered_dars"
LOGGER="wandb"

BATCH_SIZE=8 ###  ! ! ! !! !!
SEQ_LEN=2114
MAX_SHIFT=3
REVERSE_COMPLEMENT=true

LEARNING_RATE=0.00001 ### !!!!!!!
MAX_WEIGHT=100
LOSS_MULTIPLIER=1.0

FIRST_CONV_FILTERS=512
NUM_FILTERS=512

EPOCHS=100
LR_PATIENCE=3
EARLY_STOPPING_PATIENCE=10

SEED=7

PRETRAINED_MODEL_PATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/astro_human/cluster_full_universe_base_subset_filtered/checkpoints/02.keras"

#PRETRAINED_MODEL_PATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/astrocytes/astro/astro_human/cluster_full_universe_base/checkpoints/12.keras"

#PRETRAINED_MODEL_PATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/evogen-dnaseq-modeling/basal_ganglia/astro/astro_human/full_universe_base_subset_filtered_ft_2/checkpoints/07.keras"

#PRETRAINED_MODEL_PATH="/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/evogen-dnaseq-modeling/basal_ganglia/astro/astro_human/maxweight1_full_universe_base/checkpoints/17.keras"  # leave empty if unused

# === ENVIRONMENT SETUP ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CONDA_ENV

# === RUN SCRIPT ===
CMD="python /home/niklas.kempynck/nkemp/software/HMBA_Genomics/SpinalCord/crested_analysis/train_crested_model.py \
    --genome_fasta \"$GENOME_FASTA\" \
    --chr_sizes \"$CHR_SIZES\" \
    --adata_path \"$ADATA_PATH\" \
    --project_name \"$PROJECT_NAME\" \
    --run_name \"$RUN_NAME\" \
    --logger \"$LOGGER\" \
    --batch_size $BATCH_SIZE \
    --seq_len $SEQ_LEN \
    --max_stochastic_shift $MAX_SHIFT \
    --always_reverse_complement $REVERSE_COMPLEMENT \
    --learning_rate $LEARNING_RATE \
    --max_weight $MAX_WEIGHT \
    --loss_multiplier $LOSS_MULTIPLIER \
    --epochs $EPOCHS \
    --lr_patience $LR_PATIENCE \
    --early_stopping_patience $EARLY_STOPPING_PATIENCE \
    --first_conv_filters $FIRST_CONV_FILTERS \
    --num_filters $NUM_FILTERS \
    --seed $SEED"

if [ -n "$PRETRAINED_MODEL_PATH" ]; then
    CMD="$CMD --pretrained_model_path \"$PRETRAINED_MODEL_PATH\""
fi

echo "Running training script:"
echo $CMD

eval $CMD
