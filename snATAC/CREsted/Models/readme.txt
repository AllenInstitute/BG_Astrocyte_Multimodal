CREsted models trained in "Circuit-dependent specialization of human basal ganglia astrocytes"

1. Subgroup-level model
This model is the main model used in the paper predicting chromatin accessibilty for the three identified astroctyic subgroups: GM STR, GM exSTR and WM.

The model: DeepHumanAstroBG.keras
The adata file (containing the identified subgroup-level DARs and aggregrated peak heights and class names): subgroup_dars.h5ad

2. Cluster-level model
This model was trained for making cluster-level predictions within the GM STR subgroup for the following clusters: Human 14, Human 31 and Human 7 (that is h-14, h-31, and h-7).

The model: clusterlevel_DeepHumanAstroSTR.keras
The adata file (contain the identified cluster-level DARs and aggregrated peak heights and class names): cluster_dars.h5ad
