import sys
import os
import csv
import scipy
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.svm import LinearSVC,SVC
from sklearn.linear_model import SGDClassifier, LogisticRegressionCV, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, make_scorer, accuracy_score, classification_report, confusion_matrix
from sklearn.utils.multiclass import unique_labels
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit, GridSearchCV
import scanorama as scan 


def read_data(train_paths, train_datasets, test_path, test_dataset, use_raw = False):

    """ Function to read in  training and testing datasets
    
    This function reads a list of training datasets (at least one) and one testing dataset, either from raw or processed. It returns a list of training dataset adata object, and the testing adata object.
    
    parameters
    ----------
    train_paths: `list`
        list of paths where training datasets are located
    train_datasets: `list`
        a list of names of training datasets (same order as paths) 
    test_path: `string` 
        path of test dataset
    test_dataset: `string`
        name of test dataset
    use_raw: `bool` | default = False
        whether to use processed adata object or cp10k normalized counts

    returns
    -------
    `list`
        lits of training anndata objects.

    AnnData
        Anndata object with testing dataset
    AnnData
        test dataset anndata object that will not be modified, is used to append prediction column to.
    """

    adata_ignore, adata_orig = read_adata(train_paths, train_datasets, test_path, test_dataset)

    if use_raw == True:
        print('Reading files from raw')
        adata_trains, adata_pred = read_raw(train_paths, train_datasets, test_path, test_dataset)
    else:
        print('Reading files')
        #return(train_datasets)
        adata_trains, adata_pred = read_adata(train_paths, train_datasets, test_path, test_dataset)
    return adata_trains, adata_pred, adata_orig

def read_raw(train_paths, train_datasets, test_path, test_dataset):
    """ read from adata.raw and revert log1p normalization
    
    This function reads a list of training datasets (at least one) and one testing dataset and reverses log1p normalization on the raw set that includes all genes and has not been regressed out

    parameters
    ----------
    train_paths: `list`
        list of paths where training datasets are located
    train_datasets: `list`
        a list of names of training datasets (same order as paths) 
    test_path: `string` 
        path of test dataset
    test_dataset: `string`
        name of test dataset
    
    returns
    -------
    list
        A list of training dataset anndata objects

    AnnData
        An anndata object containing testing dataset
    """

    adata_trains =[]
    for i in range(len(train_datasets)):
        adata_trains.append(scv.read(os.path.join(train_paths[i], train_datasets[i])))
        adata_trains[i]= sc.AnnData(X = np.expm1(adata_trains[i].raw.X), obs = adata_trains[i].obs, var= adata_trains[i].raw.var)
    adata_pred =scv.read(os.path.join(test_path, test_dataset))
    adata_pred = sc.AnnData(X = np.expm1(adata_pred.raw.X), obs = adata_pred.obs, var= adata_pred[i].raw.var)
    
    return adata_trains, adata_pred

def read_adata(train_paths, train_datasets, test_path, test_dataset):
    """ read adata files of training and testing datasets
    
    This function reads a list of training datasets (at least one) and one testing dataset from .h5ad files and returns a list of training datasets anndata objects, and a testing anndata object.
    
    parameters
    ----------
    train_paths: `list`
        list of paths where training datasets are located
    train_datasets: `list`
        a list of names of training datasets (same order as paths) 
    test_path: `string` 
        path of test dataset
    test_dataset: `string`
        name of test dataset
    

    returns
    -------
    list
        A list of training dataset anndata objects

    AnnData
        An anndata object containing testing dataset
    """

    adata_trains =[]
    for i in range(len(train_datasets)):
        adata_trains.append(scv.read(os.path.join(train_paths[i], train_datasets[i])))
    adata_pred =scv.read(os.path.join(test_path, test_dataset))
    return adata_trains, adata_pred

def merge_data(adata_trains, adata_pred, genes_to_use = 'all', merge = 'scanorama'):

    """ read adata files of training and testing datasets
    
    This function reads a list of training datasets (at least one) and one testing dataset from .h5ad files and returns a list of training datasets anndata objects, and a testing anndata object.
    
    parameters
    ----------
    adata_trains: `list`
        list of training adata objects
    adata_pred: `list`
        testing adata object
    train_datasets: `list` 
        list of name of training datasets
    genes_to_use: `list` or `string` | default = 'all'
        if `all` nothing happens, otherwise all genes not found in the list will be removed.
    merge: `string` | default = 'scanorama'
        merges datasets using scanorama. if time is an issue, choose 'naive' for simple concatenation.
   

    returns
    -------
    list
        A merged and corrected training adata object containing chosen genes

    AnnData
        An anndata object containing corrected testinf adata object with chosen genes
    """

    if genes_to_use != 'all':
        adata_trains, adata_pred = remove_genes(adata_trains, adata_pred, genes_to_use)
    
    # merge datasets if there are multiple training_datasets
    
    if merge == 'naive':
        print('merging naively')
        adata_train = naive_merge(adata_trains, train_datasets)
    # implement scanorama and scanorama integrated
    if merge == 'scanorama':
        print('merging with scanorama')
        adata_train, adata_pred = scanorama_merge(adata_trains, adata_pred, True)
    else:
        adata_train = adata_trains[0]
    if merge == 'naive' or len(adata_trains) ==1:
        print('calculating intersection')
        adata_train, adata_pred = intersect_genes(adata_train, adata_pred)    
    return adata_train, adata_pred

def naive_merge(adata_trains):
    """ concatenates training anndata objects
    
    parameters
    ----------
    adata_trains: `list`
        list of training anndata objects
    train_datasets: `list`
        a list of names of training datasets (same order as adata_trains) 

    returns
    -------

    AnnData
        A concatenated anndata object of all training datasets.
    """
    # note: the following lines are a hack to avoid an incomprehensible error when aome values i uns are not the same
    """
    rmkeys = ['neighbors', 'pca', 'rank_genes_groups', 'celltype_dream_colors']
    for trains in adata_trains:
        # note: anndata concatenate only works if dtypes of var are the same
        trains.var.SYMBOL = trains.var.SYMBOL.astype('object')
        if trains.raw:
            trains.raw.var.SYMBOL = trains.raw.var.SYMBOL.astype('object')
        for key in rmkeys:
            trains.uns.pop(key, None)
            trains.uns.pop(key, None)
    """

    if len(adata_trains) ==1:
        return adata_trains[0]
    adata_train = adata_trains[0].concatenate(*adata_trains[1:])
    return adata_train


def scanorama_merge(adata_trains, adata_pred, keepdimensionality):
    """ corrects datasets using scanorama and merge training datasets subsequently
    
    This function reads a list of training datasets (at least one) and one testing dataset from .h5ad files and returns a merged and corrected training dataset anndata object, and a corrected testing anndata object.
    
    parameters
    ----------
    adata_trains: `list`
        list of training dataset adata objects
    adata_pred: AnnData
        testing dataset anndata object
    keepdimensionality: `bool` 
        determines if we should use all common genes or if we should reduce dimensionality to 100. False not currently implemented    
    train_datasets: `list`
        names of train datasets

    returns
    -------
    AnnData
        A concatenated and corrected anndata object of all training datasets.
    AnnData
        An anndata object containing corrected testing dataset.
    """

    adata_pred_obssave = adata_pred
    nonmerged_adata_train = naive_merge(adata_trains) # to have merged obs ready
    all_adata = naive_merge([nonmerged_adata_train, adata_pred])
    adata_trains.append(adata_pred)
    print('using scanorama rn')
    integrated, corrected = scan.correct_scanpy(adata_trains, return_dimred=True) 
    print('integrating training set')
    if len(adata_trains) != 2:
        adata_train = sc.AnnData.concatenate(*corrected[:-1])
    else:
        adata_train = corrected[0]
    adata_train.obs = nonmerged_adata_train.obs
    adata_train.var = all_adata.var
    adata_pred = sc.AnnData(corrected[-1])
    adata_pred.obs = adata_pred_obssave.obs
    adata_pred.var = all_adata.var

    adata_trains.pop()
    return adata_train, adata_pred


def remove_genes(adata_trains, adata_pred, genes_to_use):
    """ removes all genes not in gene set
    
    removes all genes from var that are not found in specified gene set list
    
    parameters
    ----------
    adata_trains: `list`
        list of training dataset adata objects
    adata_pred: AnnData
        testing dataset anndata object
    genes_to_use: `list` 
        list of genes to be removed
    

    returns
    -------
    list
        A list of training dataset anndata objects

    AnnData
        An anndata object containing testing dataset
    """

    for i in range(len(adata_trains)):
        adata_trains[i] = adata_trains[i][:,adata_trains[i].var['SYMBOL'].isin(genes_to_use)]
    adata_pred = adata_pred[:,adata_pred.var['SYMBOL'].isin(genes_to_use)]
    return adata_trains, adata_pred




def intersect_genes(adata_train, adata_pred):
    """ removes all genes not in all datasets
    
    removes all genes not contained in all datasets
    
    parameters
    ----------
    adata_train: AnnData
        training dataset anndata object
    adata_pred: AnnData
        testing dataset anndata object

    
    returns
    -------
    AnnData
        An anndata object containing training dataset with non-shared genes removed

    AnnData
        An anndata object containing testing dataset with non-shared genes removed
    """

    """Find intersecting subset of genes between sample and signature matrix"""

    # note: the following lines are a hack to avoid an incomprehensible error when aome values i uns are not the same
    rmkeys = ['neighbors', 'pca', 'rank_genes_groups', 'celltype_dream_colors']
    for key in rmkeys:
        adata_train.uns.pop(key, None)
        adata_pred.uns.pop(key, None)

    intersection = pd.concat([adata_train.var, adata_pred.var], axis = 1, join='inner') #, on = index
    intersection = intersection.index#.tolist()
    adata_train = adata_train[:,adata_train.var.index.isin(intersection)]
    adata_pred = adata_pred[:,adata_pred.var.index.isin(intersection)]
    return adata_train, adata_pred


def remove_nonshared(adata_train, adata_pred):

    """ removes all celltypes not in all datasets
    
    NOTE: requires annotated testing dataset
    
    parameters
    ----------
    adata_train: AnnData
        training dataset anndata object
    adata_pred: AnnData
        testing dataset anndata object

    
    returns
    -------
    AnnData
        An anndata object containing training dataset with non-shared celltypes removed

    AnnData
        An anndata object containing testing dataset with non-shared celltpyes removed
    """
    shared = pd.Series(list(set(adata_train.obs[celltype]) & set(adata_pred.obs[celltype])))
    adata_train_shared = adata_train[adata_train.obs[celltype].isin(shared)]
    adata_pred_shared = adata_pred[adata_pred.obs[celltype].isin(shared)]
    removed_train = adata_train.obs.shape[0] - adata_train_shared.obs.shape[0]
    removed_test = adata_pred.obs.shape[0] - adata_pred_shared.obs.shape[0]
    remove_train_percentage = removed_train/ adata_train.obs.shape[0] 
    remove_test_percentage = removed_test/ adata_pred.obs.shape[0] 
    return adata_train_shared, adata_pred_shared


def fit(adata_train, method, celltype, njobs =10):
    """ fits classifier on training dataset
    
    uses specified celltype column as label
    
    parameters
    ----------
    adata_train: AnnData
        training dataset anndata object
    method: `string`
        either 'linear' for a linar svm, 'sgd' for stochastic gradient descent, recommended when using raw, 'rbf' radial basis function, not recommended due to time
        it is also possible to choose a random forest classifier with 'random_forest', is faster than other methods but performs less well, or a multiclass logisitc regression using 'logistic_regression'. this allows you, like the random forest to specify cutoffs and have cells classifier as unknowns. Recommended.
    celltype: `string`
        column name of column to be used for classification
    njobs: int
        number of cores to use, only applies to regression or random forest classifiers

    
    returns
    -------
    sklearn.calibration.CalibratedClassifierCV or other classifier class
        trained svm classifier
    sklearn.preprocessing.StandardScaler
        a scaler fitted on the training set to be used on testing set
    """

    if scipy.sparse.issparse(adata_train.X) == True:
        train = adata_train.X.todense()

    else:
        train = adata_train.X
    
    scaler = StandardScaler().fit(train) 
    train = scaler.transform(train)

    train = pd.DataFrame(train)
    y_train = adata_train.obs[celltype].ravel()
    
    if method == 'linear':
        classifier = linear_svm(train, y_train)

    if method == 'rbf':
        classifier = rbf_svm(train, y_train)

    if method == 'sgd':
        classifier = sgd_svm(train, y_train)

    if method == 'random_forest':
        classifier = random_forest(train, y_train, njobs)

    if method == 'logistic_regression':
        classifier = logistic_regression(train, y_train, njobs) 

    if method == 'logistic_regression_ovr':
        classifier = logistic_regression_ovr(train, y_train, njobs) 

    if method == 'logistic_regression_elastic':
        classifier = logistic_regression_elastic(train, y_train, njobs)

    return classifier, scaler




def linear_svm(train, y_train):

    """ fits linear svm on training dataset
    
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframeelastic
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.calibration.CalibratedClassifierCV
        trained svm classifier
    """
    
    svm = LinearSVC()
    calib_svm = CalibratedClassifierCV(svm)
    calib_svm.fit(train, y_train)
    return calib_svm

def rbf_svm(train, y_train):

    """ fits radial basis function kernel svm on training dataset
    
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframe
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.model_selection.GridSearchCV
        trained svm classifier
    """
    C_range = np.logspace(-2, 10, 13)
    gamma_range = np.logspace(-9, 3, 13)  
    param_grid = dict(gamma=gamma_range, C=C_range)      
    cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)             
    GS_svm = GridSearchCV(SVC(), param_grid=param_grid, cv=cv, n_jobs = 20) 
    GS_svm.fit(train, y_train)
    return GS_svm

def sgd_svm(train, y_train):
    """ fits linear svm on training dataset using stochastic gradient descent
    
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframe
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.calibration.CalibratedClassifierCV
        trained svm classifier
    """

    svm = SGDClassifier(max_iter=1000, tol = 1e-3)
    svm.fit(train, y_train)
    calib_svm = CalibratedClassifierCV(svm)
    calib_svm.fit(train, y_train)
    return calib_svm


def random_forest(train, y_train, njobs):
    """ fits a random forest of a thousand esitamtors with balance class weight on training dataset. note: need at least ten nodes
    
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframe
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.ensemble.RandomForestClassifier
        trained svm classifier
    """

    clf = RandomForestClassifier(n_estimators= 1000, random_state=0, verbose = 1,class_weight = 'balanced',  n_jobs = njobs)
    clf.fit(train, y_train)
    
    return clf


def logistic_regression(train, y_train, njobs):
    """ multiclass crossvalidated logistic regression with balanced class weight. Requires ten nodes.
    
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframe
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.calibration.CalibratedClassifierCV
        trained svm classifier
    """
    clf = LogisticRegressionCV(random_state=0, verbose= 1, class_weight = 'balanced', n_jobs = njobs, multi_class= 'multinomial').fit(train, y_train) 
    return clf


def logistic_regression_ovr(train, y_train, njobs):
    """ multiclass crossvalidated logistic regression with balanced class weight. Requires ten nodes. Returns OVR probability.
    
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframe
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.calibration.CalibratedClassifierCV
        trained
    """
    clf = LogisticRegressionCV(random_state=0, verbose= 1, class_weight = 'balanced',  n_jobs = njobs, multi_class= 'ovr').fit(train, y_train)

    return clf

def logistic_regression_elastic(train, y_train, njobs):
    """ multiclass crossvalidated logistic regression with balanced class weight. Requires ten nodes. Uses elastic regularisation.
    parameters
    ----------
    train: pd.DataFrame
        adata_train.X but scaled and as a dataframe
    y_train: pd.DataFrame
            one dimensional dataframe containing class label
    
    returns
    -------
    sklearn.calibration.CalibratedClassifierCV
        trained
    """
    clf = LogisticRegressionCV(penalty = 'elasticnet', solver = 'saga', random_state=0, verbose= 1, class_weight = 'balanced',  n_jobs = njobs, multi_class= 'multinomial',  cv = 3, l1_ratios = [0.15]).fit(train, y_train)

    return clf



def adata_predict(classifier, scaler,adata_pred, adata_orig, threshold = 0):
    """ predicts on testing set using trained classifier
    
    parameters
    ----------
    test_path: `string` 
        path of test dataset
    test_dataset: `string`
        name of test dataset
    adata_orig: AnnData
        Test dataset can also be provided as an anndata object. It is different from adata_pred in that it has not undergone scanorama merging.
    classifier: sklearn object
        trained classifier of chosen type
    scaler: sklearn.preprocessing.StandardScaler
        scaler fit to training dataset
    adata_pred: AnnData
        Preprocessed anndata object containing testing dataset.
    threshold: `float` | default = 0
        value between 0 and 1, used for classifying into unknowns

    returns
    -------
    AnnData
        Original anndata object with annotation added as column
    """
    pred = predict(classifier, scaler, adata_pred, threshold)
    adata_orig.obs['auto_annot'] = pred

    return adata_orig

def predict(classifier, scaler, adata_pred, threshold = 0):

    """ predicts on testing set using trained classifier
    
    parameters
    ----------
    classifier: sklearn object
        trained classifier of chosen type
    scaler: sklearn.preprocessing.StandardScaler
        scaler fit to training dataset
    adata_pred: AnnData
        Preprocessed anndata object containing testing dataset.
    threshold: `float` | default = 0
        value between 0 and 1, used for classifying into unknowns

    returns
    -------
    np.array
        array of most likely class label for every cell
    """

    if scipy.sparse.issparse(adata_pred.X) == True:
        test = adata_pred.X.todense()
    else:
        test = adata_pred.X
    test = scaler.transform(test)
    results = classifier.predict(test)

    if (isinstance(classifier, LogisticRegressionCV) and classifier.multi_class == 'ovr') == False:
        prob = np.max(classifier.predict_proba(test), axis=1)

    else:  # in this case we are using ovr logreg
        prob = np.max(scipy.special.expit(classifier.decision_function(test)), axis=1)

    unlabeled = np.where(prob < threshold)
    results[unlabeled] = 'unknown'
    pred = pd.DataFrame(results)

    pred.to_csv("SVM_Pred_Labels_inter_jupyter.csv", index = False)
    return results 


def report(adata_pred, celltype, method, analysis_name, train_datasets, test_dataset, merge, use_raw = False, genes_to_use = 'all', remove_nonshared = False, clustering = 'leiden', asymmetric_matrix = True):
    """ reports basic metrics, produces confusion matrices and plots umap of prediction

    Writes out a csv file containing all accuracy and f1 scores.
    Writes normalized and absolute confusion matrices, as well as umap prediction comparisons to ./figures.
    
    parameters
    ----------
    adata_pred: AnnData
        original adata object with 'auto_annot' column
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
    0
    """

    # calculate umaps for plot
    if 'X_umap' not in adata_pred.obsm:
        sc.tl.umap(adata)


    #get acc
    acc = accuracy_score(adata_pred.obs[celltype], adata_pred.obs['auto_annot'])

    #get f1
    f1 = f1_score(adata_pred.obs[celltype], adata_pred.obs['auto_annot'], labels = adata_pred.obs[celltype], average = 'macro')

    #get report
        
    report = classification_report(adata_pred.obs[celltype], adata_pred.obs['auto_annot'], output_dict=True)
    sklearn_report = pd.DataFrame(report).transpose()
    
    #csv file with important metrics
    with open('auto_annot_report_'  + analysis_name +'.csv', mode='w') as report:
        report_writer = csv.writer(report, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        report_writer.writerow(['train_dataset = ', train_datasets, 'test_datset = ', test_dataset])
        report_writer.writerow(['celltype = ', celltype, 'method = ', method])
        report_writer.writerow(['remove_nonshared = ', remove_nonshared, 'merge = ', merge, 'use_raw = ', use_raw, 'genes_to_use = ', genes_to_use])

        report_writer.writerow(['accuracy=', acc, 'f1=', f1])

        report_writer.writerow(['classification report'])
        sklearn_report.to_csv(report, header = True)
        
        
        # make umap
        sc.settings.set_figure_params(dpi=240)

        sc.pl.umap(adata_pred, color=[celltype, 'auto_annot', clustering], legend_loc='on data',legend_fontsize=7, save= '.ondata_'+ analysis_name + '.png')
        sc.pl.umap(adata_pred, color=[celltype, 'auto_annot', clustering],legend_fontsize=7, wspace = 1.5,  save = '.' + analysis_name + '.png')
        sc.settings.set_figure_params(dpi=60)

        # make conf matrices (4)
        class_names =  np.unique(np.concatenate((adata_pred.obs[celltype], adata_pred.obs['auto_annot'])))
        np.set_printoptions(precision=2)
        # Plot non-normalized confusion matrix
        plot_confusion_matrix(adata_pred.obs[celltype], adata_pred.obs['auto_annot'], classes=class_names, celltype=celltype,  title='Confusion matrix, without normalization',numbers = False, adata_predicted = adata_pred, asymmetric_matrix = asymmetric_matrix)
        plt.savefig(os.path.join('./figures/SVM_confusion_matrix_' + analysis_name +'_' + celltype +  '.svg'))

        # Plot normalized confusion matrix with numbers
        plot_confusion_matrix(adata_pred.obs[celltype], adata_pred.obs['auto_annot'], classes=class_names,celltype=celltype,  normalize=True, title='Normalized confusion matrix', numbers = False, adata_predicted = adata_pred, asymmetric_matrix = asymmetric_matrix)
        plt.savefig(os.path.join('./figures/SVM_confusion_matrix_norm_' + analysis_name +'_' + celltype +  '.svg'))

    
def plot_confusion_matrix(y_true, y_pred, classes, celltype,
                          normalize=False,
                          title=None, numbers =False,
                          cmap=plt.cm.Blues, adata_predicted= None, asymmetric_matrix = True): 
    """ plots confusion matrices

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
    matplotlib.use('Agg')
    
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion matrix, without normalization'

    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    # Only use the labels that appear in the data
    #classes = classes[unique_labels(y_true, y_pred)]
    if asymmetric_matrix == True:
        class_names =  np.unique(np.concatenate((adata_predicted.obs[celltype], adata_predicted.obs['auto_annot'])))
        class_names_orig = np.unique(adata_predicted.obs[celltype])
        class_names_pred = np.unique(adata_predicted.obs['auto_annot'])
        test_celltypes_ind = np.searchsorted(class_names, class_names_orig)
        train_celltypes_ind = np.searchsorted(class_names, class_names_pred)
        cm=cm[test_celltypes_ind,:][:,train_celltypes_ind]
    
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    fig, ax = plt.subplots(figsize=(15,15))
    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    ax.figure.colorbar(im, ax=ax, shrink = 0.8)
    # We want to show all ticks...
    if asymmetric_matrix == True:
        ax.set(xticks=np.arange(cm.shape[1]),
               yticks=np.arange(cm.shape[0]),
               # ... and label them with the respective list entries
               xticklabels=class_names_pred, yticklabels=class_names_orig,
               title=title,
               ylabel='True label',
               xlabel='Predicted label')
    else:
        ax.set(xticks=np.arange(cm.shape[1]),
               yticks=np.arange(cm.shape[0]),
               # ... and label them with the respective list entries
               xticklabels=classes, yticklabels=classes,
               title=title,
               ylabel='True label',
               xlabel='Predicted label')
        
    ax.grid(False)
    #ax.tick_params(axis='both', which='major', labelsize=10)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    if numbers == True:
        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, format(cm[i, j], fmt),
                        ha="center", va="center",
                        color="white" if cm[i, j] > thresh else "black")
    #fig.tight_layout()
    return ax

def adata_pred_prob(classifier, scaler,adata_pred, adata_orig, threshold = 0):
    """ predicts on testing set using trained classifier and returns class probability for every cell and every class
    
    parameters
    ----------
    test_path: `string` 
        path of test dataset
    test_dataset: `string`
        name of test dataset
    classifier: sklearn object
        trained classifier of chosen type
    scaler: sklearn.preprocessing.StandardScaler
        scaler fit to training dataset
    adata_pred: AnnData
        Preprocessed anndata object containing testing dataset.
    threshold: `float` | default = 0
        value between 0 and 1, used for classifying into unknowns

    returns
    -------
    AnnData
        Original anndata object with annotation added as column and a column for each training celltype probability
    """
    adata_prob = adata_orig.copy()
    pred, prob = predict_proba(classifier, scaler, adata_pred, threshold)

    adata_prob.obs['auto_annot'] = pred
    adata_prob.obs[classifier.classes_] = pd.DataFrame(prob,  index=adata_prob.obs.index)
    return adata_prob

def predict_proba(classifier, scaler, adata_pred, threshold = 0):

    """ predicts on testing set using trained classifier and returns probabilities
    
    parameters
    ----------
    classifier: sklearn object
        trained classifier of chosen type
    scaler: sklearn.preprocessing.StandardScaler
        scaler fit to training dataset
    adata_pred: AnnData
        Preprocessed anndata object containing testing dataset.
    threshold: `float` | default = 0
        value between 0 and 1, used for classifying into unknowns

    returns
    -------
    np.array
        array of most likely class label for every cell
    np.array
        array of all probabilities of all cell types
    """

    if scipy.sparse.issparse(adata_pred.X) == True:
        test = adata_pred.X.todense()
    else:
        test = adata_pred.X
    test = scaler.transform(test)

    results = classifier.predict(test)

    if (isinstance(classifier, LogisticRegressionCV) and classifier.multi_class == 'ovr') == False:
        probmax = np.max(classifier.predict_proba(test), axis=1)
        prob = classifier.predict_proba(test)
        # (pd.DataFrame(prob)) #test if this is needed

    else:  # in this case we are using ovr logreg
        probmax = np.max(scipy.special.expit(classifier.decision_function(test)), axis=1)
        prob = (scipy.special.expit(classifier.decision_function(test)))
        # (pd.DataFrame(prob))

    unlabeled = np.where(probmax < threshold)
    results[unlabeled] = 'unknown'
    pred = pd.DataFrame(results)

    pred.to_csv("SVM_Pred_Labels_inter_jupyter.csv", index=False)

    return results, prob 
