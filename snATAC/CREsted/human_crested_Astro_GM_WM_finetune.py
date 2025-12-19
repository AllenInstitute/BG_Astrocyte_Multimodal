import os, re, sys, glob
import anndata as ad
import pandas as pd
import keras
import crested
import pysam
from tqdm import tqdm

import tensorflow as tf

gpus = tf.config.list_physical_devices('GPU')
if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
        print("Memory growth enabled for GPUs")
    except RuntimeError as e:
        print("Failed to set memory growth:", e)

## --------------------------
## Load model
species_name = "human"
work_dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling"

## -------------------------
## Setup genome
reference = os.path.join(work_dir, f"genomes/{species_name}")
if species_name == "macaque":
    reference = os.path.join(reference, "ncbi")

## Chromosome sizes
if os.path.exists(os.path.abspath(os.path.join(reference, f"{species_name}.chrom.sizes"))):
    chr_sizes = os.path.abspath(os.path.join(reference, f"{species_name}.chrom.sizes"))
else:
    chr_sizes = None

## Fasta and GTF paths
fasta_path = os.path.abspath(os.path.join(reference, f"fasta/genome.fa"))
gtf_path = os.path.abspath(os.path.join(reference, f"genes/genes.gtf.gz"))

##
genome = crested.Genome(
    fasta=fasta_path,
    annotation=gtf_path,
)
crested.register_genome(
    genome
)

##
best_model = "/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/evogen-dnaseq-modeling/basal_ganglia/CrestedModelAPI--dilated-cnn-human-basalganglia-astro-gm-wm--25-06-08-02-18-44/checkpoints/20.keras"

## Initalize model
model_architecture = keras.models.load_model(best_model, compile=False)

## ---------------------------
## Gather data model was trained on, really only needed to determine class (celltype) names
adata = ad.read_h5ad(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata", f"{species_name}_basalganglia_hmba_Astro_GM_WM_pre-print_crested.h5ad"))

## All regions with a Gini index 1 std above the mean across all regions will be kept
crested.pp.filter_regions_on_specificity(
    adata, gini_std_threshold=1.0
) 
adata.write_h5ad(os.path.join(work_dir, "basal-ganglia", "data", species_name, "crested_adata", f"{species_name}_basalganglia_hmba_Astro_GM_WM_pre-print_crested_finetune.h5ad"))

## 
datamodule = crested.tl.data.AnnDataModule(
    adata,
    genome=genome,
    batch_size=64,  ## Recommended to go for a smaller batch size than in the pretrained model
    max_stochastic_shift=3,
    always_reverse_complement=True,
    chromsizes_file=chr_sizes,
)

## Use the same config you used for the pretrained model. EXCEPT THE LEARNING RATE, 
## make sure that is lower than it was on the epoch you select the model from.
optimizer = keras.optimizers.Adam(learning_rate=1e-5)  ## Lower LR!
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100)
metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]

alternative_config = crested.tl.TaskConfig(optimizer, loss, metrics)
print(alternative_config)

## setup the trainer
os.chdir("/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/OCTO/aibs-octo-dnaseq-modeling/basal-ganglia/models/crested/evogen-dnaseq-modeling/basal_ganglia/CrestedModelAPI--dilated-cnn-human-basalganglia-astro-gm-wm--25-06-08-02-18-44/") ## Kind of annoying that CREsted operates in a specific directory.
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=alternative_config,
    project_name=species_name,  # change to your liking
    run_name="checkpoints_finetune",  # change to your liking
    logger="tensorboard",  # or 'wandb', 'tensorboard'
)

## Train!
trainer.fit(epochs=60,
            model_checkpointing_best_only=True)
