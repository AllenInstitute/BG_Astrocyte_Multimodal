# import statements
import pandas as pd
import matplotlib.pyplot as plt

def plot_confusion_matrix(
    obs_df, 
    annot_col, 
    mapped_col, 
    cmap='inferno',
    annot_values_on_heatmap=True,
    fig_size=(4,4),
    vlims=None,
):
    """
    Plot raw and normalized confusion matrices for two categorical columns.
    
    Parameters
    ----------
    obs_df : pd.DataFrame
        The .obs dataframe containing the columns to compare.
    annot_col : str
        Column name to use as true labels (plotted on columns).
    mapped_col : str
        Column name to use as mapped labels (plotted on rows).
    cmap : str, optional
        Colormap for the heatmaps. Default is 'magma_r'.
    """
    import seaborn as sns
    from sklearn.metrics import confusion_matrix

    true_labels = obs_df[annot_col]
    mapped_labels = obs_df[mapped_col]

    # # all_labels = sorted(set(true_labels.unique()) | set(mapped_labels.unique()))
    # all_labels = list(pd.Categorical(true_labels).categories)

    # cm = confusion_matrix(true_labels, mapped_labels, labels=all_labels)

    if hasattr(true_labels, 'cat'):
        all_labels = list(true_labels.cat.categories)
    else:
        all_labels = list(pd.Categorical(true_labels).categories)
    
    true_labels = true_labels.astype(str)
    mapped_labels = mapped_labels.astype(str)
    all_labels = [str(label) for label in all_labels]

    cm = confusion_matrix(true_labels, mapped_labels, labels=all_labels)

    cm_df = pd.DataFrame(
        cm.T,
        index=all_labels,
        columns=all_labels
    )
    cm_df.index.name = f'Mapped ({mapped_col})'
    cm_df.columns.name = f'Annotated ({annot_col})'

    # cm_norm = cm_df.div(cm_df.sum(axis=0), axis=1)
    cm_norm = cm_df.div(cm_df.sum(axis=0), axis=1).fillna(0)

    # n_labels = len(all_labels)

    fig, ax = plt.subplots(1, 1, figsize=fig_size)

    if vlims is not None:
        vmin, vmax = vlims
    else:
        vmin, vmax = None, None

    # Create a custom annotation matrix for non-zero values
    if annot_values_on_heatmap:
        annot_matrix = cm_norm.apply(lambda row: row.map(lambda x: f"{x:.2f}" if x >= 0.005 else ""))   
        fmt = ''  # No need for additional formatting since values are already formatted in the matrix
    else:
        annot_matrix = None
        fmt = '.2f'  # Default formatting for any annotations if they were to be added

    sns.heatmap(
        cm_norm,
        ax=ax,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        annot=annot_matrix,
        fmt=fmt,
        square=True,
        cbar_kws={'label': 'Fraction Mapped'},
        xticklabels=True,
        yticklabels=True
    )
    ax.set_title(f'Confusion Matrix: {annot_col} vs {mapped_col}\n(normalized by annotated label)')
    ax.tick_params(axis='x', rotation=90)
    ax.tick_params(axis='y', rotation=0)

    plt.tight_layout()
    plt.show()

    return fig