import csv
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
    make_scorer,
)

from ..pl._riverplot import riverplot_2categories




def report(
    adata_pred,
    celltype,
    method,
    analysis_name,
    train_datasets,
    test_dataset,
    merge,
    name_prediction="auto_annot",
    name_report="auto_annot",
    use_raw=False,
    genes_to_use="all",
    remove_nonshared=False,
    clustering="leiden",
    asymmetric_matrix=True,
):
    """reports basic metrics, produces confusion matrices and plots umap of prediction

    Writes out a csv file containing all accuracy and f1 scores.
    Writes normalized and absolute confusion matrices, as well as umap prediction comparisons to ./figures.

    parameters
    ----------
    adata_pred: AnnData
        original adata object with name_prediction column
    celltype: `str`
        celltype column on which the prediction was performed
    method: `str`
        method that was used for prediction.
    analysis_name: `str`
        name of the analyis, used for writing files
    train_datasets: `list`
        list of used training datasets
    test_dataset: `str`
        name of test dataset
    merge: `str`
        what merging was performed
    name_prediction : "auto_annot"| default = "auto_annot"
        observation name containing the prediction to compare with.
    name_report : "auto_annot"| default = "auto_annot"
        prefix of the report
    use_raw: `bool`  | default = False
        if anndata.raw was used
    genes_to_use: `list` or `string` | default = 'all'
        what geneset wsa used
    remove_nonshared: `bool`|default = False
    clustering: `str` | default = leiden
        clustering that was used in original analysis of testing set, needed for umap plotting
    asymmetric_matrix: `bool` | default = True
        if False returns square confusion matrix, if True it only shows possible combinations

    returns
    -------
    Figure
        A matplotlib figure element containing the riveplot generated for interactive display.
        
    """

    # calculate umaps for plot
    if "X_umap" not in adata_pred.obsm:
        sc.tl.umap(adata_pred)

    if name_prediction not in adata_pred.obs.keys():
        sys.exit(
            name_prediction
            + " label not found in the predicted dataset (should be in obs)"
        )
    # get acc
    acc = accuracy_score(adata_pred.obs[celltype], adata_pred.obs[name_prediction])

    # get f1
    f1 = f1_score(
        adata_pred.obs[celltype],
        adata_pred.obs[name_prediction],
        labels=adata_pred.obs[celltype],
        average="macro",
    )

    # get report

    report = classification_report(
        adata_pred.obs[celltype], adata_pred.obs[name_prediction], output_dict=True
    )
    sklearn_report = pd.DataFrame(report).transpose()

    # csv file with important metrics
    with open(name_report + "_report_" + analysis_name + ".csv", mode="w") as report:
        report_writer = csv.writer(
            report, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )

        report_writer.writerow(
            ["train_dataset = ", train_datasets, "test_datset = ", test_dataset]
        )
        report_writer.writerow(["celltype = ", celltype, "method = ", method])
        report_writer.writerow(
            [
                "remove_nonshared = ",
                remove_nonshared,
                "merge = ",
                merge,
                "use_raw = ",
                use_raw,
                "genes_to_use = ",
                genes_to_use,
            ]
        )

        report_writer.writerow(["accuracy=", acc, "f1=", f1])

        report_writer.writerow(["classification report"])
        sklearn_report.to_csv(report, header=True)

    # make umap
    sc.settings.set_figure_params(dpi=240)

    sc.pl.umap(
        adata_pred,
        color=[celltype, name_prediction, clustering],
        legend_loc="on data",
        legend_fontsize=7,
        save=".ondata_" + analysis_name + ".png",
    )
    sc.pl.umap(
        adata_pred,
        color=[celltype, name_prediction, clustering],
        legend_fontsize=7,
        wspace=1.5,
        save="." + analysis_name + ".png",
    )
    sc.settings.set_figure_params(dpi=60)

    # make conf matrices (4)
    class_names = np.unique(
        np.concatenate((adata_pred.obs[celltype], adata_pred.obs[name_prediction]))
    )
    np.set_printoptions(precision=2)
    # Plot non-normalized confusion matrix
    plot_confusion_matrix(
        adata_pred.obs[celltype],
        adata_pred.obs[name_prediction],
        classes=class_names,
        celltype=celltype,
        name_prediction=name_prediction, 
        title="Confusion matrix, without normalization",
        numbers=False,
        adata_predicted=adata_pred,
        asymmetric_matrix=asymmetric_matrix,
    )
    plt.savefig(
        os.path.join(
             "./figures/" + method + "_confusion_matrix_"
            + analysis_name
            + "_"
            + celltype
            + ".svg"
        )
    )

    # Plot normalized confusion matrix with numbers
    plot_confusion_matrix(
        adata_pred.obs[celltype],
        adata_pred.obs[name_prediction],
        classes=class_names,
        celltype=celltype,
        name_prediction=name_prediction, 
        normalize=True,
        title="Normalized confusion matrix",
        numbers=False,
        adata_predicted=adata_pred,
        asymmetric_matrix=asymmetric_matrix,
    )
    plt.savefig(
        os.path.join(
            "./figures/" + method + "_confusion_matrix_norm_"
            + analysis_name
            + "_"
            + celltype
            + ".svg"
        )
    )

    # plot basic riverplot
    fig = riverplot_2categories( adata=adata_pred, categories=[celltype, name_prediction])
    fig.show()
    return fig



def plot_confusion_matrix(
    y_true,
    y_pred,
    classes,
    celltype,
    name_prediction="auto_annot",
    normalize=False,
    title=None,
    numbers=False,
    cmap=plt.cm.Blues,
    adata_predicted=None,
    asymmetric_matrix=True,
):
    """plots confusion matrices

    returns a matplotlib confusion matrix

    parameters
    ----------
    y_true: pandas.core.series.Series

        ordered series of all true labels
    y_pred: pandas.core.series.Series
        ordered series of all predicted celltypes
    classes: numpy.ndarray
        union of true and predictable celltypes
    celltype: `str`
        celltype column on which the prediction was performed
    name_prediction : "auto_annot"| default = "auto_annot"
        observation name containing the prediction to compare with.
    normalize: `bool` | default = False
        whether to return absolute values or to value all celltypes equally
    title: `str` | default = None
        title to be given to confusion matrix figure in file.
    numbers: `bool`| default = False
        should the numbers be displayed in the plot. Note: is illegible in larger plots
    cmap: matplotlib.cm | default = plt.cm.Blues
        colour to be used for plotting
    asymmetric_matrix: `bool` | default = True
        if False returns square confusion matrix, if True it only shows possible combinations

    returns
    -------
    matplotlib.pyplot.plot
        plot of confusion matrix
    """
    matplotlib.use("Agg")

    if not title:
        if normalize:
            title = "Normalized confusion matrix"
        else:
            title = "Confusion matrix, without normalization"

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    # classes = classes[unique_labels(y_true, y_pred)]
    if asymmetric_matrix == True:
        class_names = np.unique(
            np.concatenate(
                (adata_predicted.obs[celltype], adata_predicted.obs[name_prediction])
            )
        )
        class_names_orig = np.unique(adata_predicted.obs[celltype])
        class_names_pred = np.unique(adata_predicted.obs[name_prediction])
        test_celltypes_ind = np.searchsorted(class_names, class_names_orig)
        train_celltypes_ind = np.searchsorted(class_names, class_names_pred)
        cm = cm[test_celltypes_ind, :][:, train_celltypes_ind]

    if normalize:
        cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print("Confusion matrix, without normalization")

    fig, ax = plt.subplots(figsize=(15, 15))
    im = ax.imshow(cm, interpolation="nearest", cmap=cmap)
    ax.figure.colorbar(im, ax=ax, shrink=0.8)
    # We want to show all ticks...
    if asymmetric_matrix == True:
        ax.set(
            xticks=np.arange(cm.shape[1]),
            yticks=np.arange(cm.shape[0]),
            # ... and label them with the respective list entries
            xticklabels=class_names_pred,
            yticklabels=class_names_orig,
            title=title,
            ylabel="True label",
            xlabel="Predicted label",
        )
    else:
        ax.set(
            xticks=np.arange(cm.shape[1]),
            yticks=np.arange(cm.shape[0]),
            # ... and label them with the respective list entries
            xticklabels=classes,
            yticklabels=classes,
            title=title,
            ylabel="True label",
            xlabel="Predicted label",
        )

    ax.grid(False)
    # ax.tick_params(axis='both', which='major', labelsize=10)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    if numbers == True:
        fmt = ".2f" if normalize else "d"
        thresh = cm.max() / 2.0
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(
                    j,
                    i,
                    format(cm[i, j], fmt),
                    ha="center",
                    va="center",
                    color="white" if cm[i, j] > thresh else "black",
                )
    # fig.tight_layout()
    return ax
