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
    adjusted_mutual_info_score,
    adjusted_rand_score,
    silhouette_score,
    pair_confusion_matrix
)

from ..pl._riverplot import riverplot_2categories




def report(
    adata_pred,
    celltype,
    method,
    analysis_name,
    train_datasets=[],
    test_dataset="",
    merge="",
    name_prediction="auto_annot",
    name_report="auto_annot",
    use_raw=False,
    genes_to_use="",
    remove_nonshared=False,
    clustering="leiden",
    asymmetric_matrix=True,
    results_folder="./",
    delimiter=",",
    verbose=False
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
    results_folder: `str` | default = './'
        output directory. A figures folder will be generated within it.
    delimiter: `str` | default = ','
        separator between fields in the csv/txt report file
    verbose: `bool` | default = False
        print verbose messages to standard out

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

    if verbose:
        print('acc: ' + str(round(acc,2)))
        print('f1: ' + str(round(f1,2)))

    # get report
    report = classification_report(
        adata_pred.obs[celltype], adata_pred.obs[name_prediction], output_dict=True
    )
    sklearn_report = round(pd.DataFrame(report).transpose(), 2)

    # get clustering scores
    ami = adjusted_mutual_info_score(adata_pred.obs[celltype], adata_pred.obs[name_prediction])
    ari = adjusted_rand_score(adata_pred.obs[celltype], adata_pred.obs[name_prediction])
    silhouette_celltype = silhouette_score(adata_pred.obsm['X_umap'], adata_pred.obs.get(celltype))
    silhouette_pred = silhouette_score(adata_pred.obsm['X_umap'], adata_pred.obs.get(name_prediction))
    pair_conf_m = pair_confusion_matrix(adata_pred.obs[celltype], adata_pred.obs[name_prediction])

    if verbose:
        print('ami: ' + str(round(ami,2)))
        print('ari: ' + str(round(ari,2)))
        print('silhouette ' + celltype + ': ' + str(round(silhouette_celltype,2)))
        print('silhouette ' + name_prediction + ': ' + str(str(round(silhouette_pred,2))))
        print('pair confusion matrix: ' + str(pd.DataFrame(pair_conf_m)))

    # csv file with important metrics
    file_ending = ".txt"
    if delimiter==",":
        file_ending = ".csv"
    with open(os.path.join(results_folder, name_report + "_report_" + analysis_name + file_ending), mode="w") as report:
        report_writer = csv.writer(
            report, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
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

        report_writer.writerow(["accuracy=", round(acc,2), "f1=", round(f1,2)])

        report_writer.writerow(["clustering report"])
        report_writer.writerow(["ari=", round(ari,2), "ami=", round(ami,2)])
        report_writer.writerow(["silhouette_celltype=", round(silhouette_celltype,2), "silhouette_pred=", round(silhouette_pred,2)])
        
        report_writer.writerow(["pair confusion matrix"])
        pd.DataFrame(pair_conf_m).to_csv(report, header=True, sep=delimiter)

        report_writer.writerow(["classification report"])
        sklearn_report.to_csv(report, header=True, sep=delimiter)

    # make umap
    sc.settings.set_figure_params(dpi=120)

    sc.pl.umap(
        adata_pred,
        color=[celltype, name_prediction, clustering],
        legend_loc="on data",
        legend_fontsize=7,
        frameon=False,
        save=".ondata_" + analysis_name + ".png",
    )
    for col in [celltype, name_prediction, clustering]:
        sc.pl.umap(
            adata_pred,
            color=col,
            wspace=1.5,
            frameon=False,
            save="." + analysis_name + "_" + col + ".png",
        )
    
    sc.settings.set_figure_params(dpi=60)
    
    os.makedirs(os.path.join(results_folder, "figures"), exist_ok=True)


    # plot basic riverplot
    riverplot = riverplot_2categories( adata=adata_pred, categories=[celltype, name_prediction])
    riverplot.show()
    riverplot.write_image(
        os.path.join(results_folder, "figures", method + "_riverplot_"
            + analysis_name
            + "_"
            + celltype
            + "_"
            + name_prediction
            + ".svg"
        )
    )
        
    # make conf matrices (4)
    class_names = np.unique(
        np.concatenate((adata_pred.obs[celltype], adata_pred.obs[name_prediction]))
    )
    np.set_printoptions(precision=2)
    # Plot non-normalized confusion matrix
    fig = plot_confusion_matrix(
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
    fig.show()
    fig.savefig(
        os.path.join(results_folder, "figures", method + "_confusion_matrix_"
            + analysis_name
            + "_"
            + celltype
            + ".svg"
        )
    )

    # Plot normalized confusion matrix with numbers
    fig = plot_confusion_matrix(
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
    fig.show()
    fig.savefig(
        os.path.join(results_folder, "figures", method + "_confusion_matrix_norm_"
            + analysis_name
            + "_"
            + celltype
            + ".svg"
        )
    )

    return riverplot



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
    #    print("Normalized confusion matrix")
    #else:
    #    print("Confusion matrix, without normalization")

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
    
    return fig
