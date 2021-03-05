from pandas import DataFrame, melt
from natsort import natsorted
from .._helper import subset_adata as _subset_adata

def count_occurrence(adata,
                    count_variable = 'celltype',
                    add_percentage = False):
    """ Generate dataframe containing the label counts/percentages of a specific column in adata.obs

    This function counts the occurrence of each label within the specified column (count_variable)
    of an AnnData object and outputs the results to a pandas DataFrame. If add_percentage is true
    it also calculates the occurrence of each label as a percentage. Note, percentages have been 
    rounded to the second decimal place after the comma.

    One of the most common use-cases for this function will be to count the occurrence of specific
    celltypes within the dataset.

    parameters
    ----------
    adata: AnnData
      the AnnData object 
    count_variable: `str` | default = 'celltype'
      string identifying the column in which the unique labels should be counted
    add_percentage: `bool` | default = False
      boolian indicator if the occurrence of each label as a percentage should be added to
      the dataframe

    returns
    -------
    pandas.DataFrame
        Dataframe containing the counts of each label, if add_percentage = True, then the DataFrame
        also contains the occurrence of each label as a percentage

    """
    #get counts for specified column
    counts = adata.obs.get(count_variable).value_counts()
    total_counts = sum(counts.tolist())    
    #initialize a dataframe to store information in
    data = DataFrame(index=counts.index.tolist(),
                     data={'Counts': counts.tolist()}  ) 
    if add_percentage:  ##  calculate percentages
        percentages = []
        for i in range(len(counts)):
            count = data['Counts'][i]
            percentages.append(round(count/total_counts*100, 2))
        #add percentages to dataframe
        data['Percentage'] = percentages    
    return(data)

def count_occurrence_subset(adata,
                           subset_variable,
                           count_variable='celltype',
                           return_percentage=False):
    """ count occurrence of a label in adata.obs after subseting adata object
    
    This function subsets the supplied AnnData object into datasubsets according to the lables contained
    in the subset_variable. For each of the datasubsets, it counts the occurrence of each label 
    within the specified column (count_variable) of the adata subset. 

    If return percentage = True then the percentage of each label occurrence within a subset is returned.
    
    parameters
    ----------
    adata: AnnData
      the AnnData object 
    subset_variable: `str`
      string identifying the column in adata.obs along which the data should be subsetted
    count_variable: `str` | default = 'celltype'
      string identifying the column in which the unique labels should be counted
    add_percentage: `bool` | default = False
      boolian indicator if the occurrence of each label as a percentage should be added to
      the dataframe

    returns
    -------
    pandas.DataFrame
        Dataframe containing the counts of each label, if add_percentage = True, then the DataFrame
        contains the occurrence of each label as a percentage within the datasubset

    """
    #get subsets
    subsets = natsorted(adata.obs.get(subset_variable).value_counts().index.tolist())

    #initialize dictionaries to store data in
    dic_counts = {}
    dic_percentages = {}
    for subset in subsets:
        #get datasubset
        data = _subset_adata(adata, filter_criteria=adata.obs.get(subset_variable) == subset, raw=False)
        all_counts = count_occurrence(data, count_variable=count_variable, add_percentage=True)

        counts = all_counts['Counts'].to_frame()
        counts.rename(columns = {'Counts':subset}, inplace=True)
        percentage = all_counts['Percentage'].to_frame()
        percentage.rename(columns = {'Percentage':subset}, inplace=True)
        
        dic_counts.update({subset:counts})
        dic_percentages.update({subset:percentage})

    #merge into one
    counts = dic_counts.get(subsets[0])
    percentages = dic_percentages.get(subsets[0])

    for subset in subsets[1:]:
      data_counts = dic_counts.get(subset)
      data_percentages = dic_percentages.get(subset)
      counts = counts.merge(data_counts, how = 'outer', left_index = True, right_index=True)
      percentages = percentages.merge(data_percentages, how = 'outer', left_index = True, right_index=True)
    
    #remove NaNs and replace with 0
    counts = counts.fillna(0)
    percentages = percentages.fillna(0)
    
    if return_percentage:
      return(percentages)
    else:
      return(counts)

def count_occurrence_subset_conditions(adata,
                                      subset_variable,
                                      condition_identifier,
                                      count_variable ='celltype',
                                      return_percentage = False):
    """count occurrence of a label for each condition in adata.obs after subseting adata object
    
    This function subsets the supplied AnnData object into datasubsets according to the lables contained
    in the subset_variable. For each of the datasubsets, it counts the occurrence of each label 
    within the specified column (count_variable) of the adata subset for each condition identified in the column condition_identifier. 

    If return percentage = True then the percentage of each label occurrence within a subset is returned.
    
    parameters
    ----------
    adata: AnnData
      the AnnData object 
    subset_variable: `str`
      string identifying the column in adata.obs along which the data should be subsetted
    count_variable: `str` | default = 'celltype'
      string identifying the column in which the unique labels should be counted
    condition_identifier: `str`
      string identifying the coloumn in which the conditions are annotated
    add_percentage: `bool` | default = False
      boolian indicator if the occurrence of each label as a percentage should be added to
      the dataframe

    returns
    -------
    pandas.DataFrame
        Dataframe containing the counts of each label, if add_percentage = True, then the DataFrame
        contains the occurrence of each label as a percentage within the datasubset

    """
    #get subsets
    subsets = adata.obs.get(subset_variable).value_counts().index.tolist()

    #get conditions
    conditions = adata.obs.get(condition_identifier).value_counts().index.tolist()
    
    #initialize dictionaries to store data in
    dic_counts={}
    dic_percentages = {}

    #generate dictionaries containing the dataframes for each subset
    for subset in subsets:
        #get datasubset
        data = _subset_adata(adata, filter_criteria = adata.obs.get(subset_variable) == subset, raw = False)

        if return_percentage:
            #generate dataframe containing percentages
            percentages = count_occurrence_subset(data, count_variable = count_variable, subset_variable = condition_identifier, return_percentage=True)
            column_names = percentages.columns.tolist()
            new_column_names = ['Percentage '+ subset + ' '+ column_name for column_name in column_names]
            percentages.columns = new_column_names

            #save to dictionary
            dic_percentages.update({subset:percentages})
        else:
        #generate dataframe containing counts
            counts = count_occurrence_subset(data, count_variable = count_variable, subset_variable = condition_identifier)
            column_names = counts.columns.tolist()
            new_column_names = ['Count '+ subset + ' '+ column_name for column_name in column_names]
            counts.columns = new_column_names

            #save to dictionary
            dic_counts.update({subset:counts})

    if return_percentage:
        #merge dictionaries into one dataframe
        percentages = dic_percentages.get(subsets[0])

        for subset in subsets[1:]:
            data_percentages = dic_percentages.get(subset)
            percentages = percentages.merge(data_percentages, how = 'outer', left_index = True, right_index=True)

        #remove NaNs and replace with 0
        percentages = percentages.fillna(0)

        #return generated dataframe
        return(percentages)
    
    else:
        counts = dic_counts.get(subsets[0])
    
        for subset in subsets[1:]:
          data_counts = dic_counts.get(subset)
          counts = counts.merge(data_counts, how = 'outer', left_index = True, right_index=True)
    
        #remove NaNs and replace with 0
        counts = counts.fillna(0)

        #return generated dataframe
        return(counts)
