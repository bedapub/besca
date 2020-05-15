#this module contains functions for cell type annotation based on signatures in python using scanpy
#import using the python version 1.3.2 at least is a prerequisit! needs to be checked in all functions

from scanpy.api import AnnData
from pandas import DataFrame, read_csv
from scipy.stats import mannwhitneyu
from numpy import log, where, array


def getset(df, signame_complete, threshold):
    """ Handles missing signatures aux function for make_anno
    Based on a dataframe of p-values, a signature name and a cutoff check if sign is present
    parameters
    ----------
    df: panda.DataFrame
      a dataframe of p-values per signature and cell cluster
    signame_complete: str
      signature name
    threshold: numpy.float64
      cutoff used for cluster attributions
    returns
    -------
    set
        Set of clusters passing the threshold value
    """

    if(bool(df.index.isin([signame_complete]).any())):
        return(set(df.columns[df.loc[signame_complete, :] > threshold]))
    else:
        return(set())


def gene_fr(f, Tcells):
    """ Returns the fraction of cells expressing a gene in distinct clusters
    Takes as an input a dataframe with fractions per clusters and a gene name
    Returns the fraction genes per cluster
    parameters
    ----------
    f: panda.DataFrame
      a dataframe of fraction genes per cell as outputted by besca
    Tcells: 'str'
      gene name

    returns
    -------
    panda.Series
        a Series of fractions per cluster and gene
    """
    k = f.loc[f['Description'].isin( Tcells),].copy()
    k = k.iloc[:,0:len(k.columns)]
    return(k.mean(axis=0, skipna=True))


def score_mw(f, mymarkers):
    """ Score Clusters based on a set of immune signatures to generate a df of pvals
    Takes as an input a dataframe with fractions per clusters and a dictionary of signatures
    Performs a Mann-Whitney test per each column and Signature, returns -10logpValues
    parameters
    ----------
    f: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    mymarkers: Dictionary
      a dictionary of signatures

    returns
    -------
    panda.DataFrame
        a dataframe of -10logpValues per cluster and signature
    """
    ids = set(f.iloc[:,2:len(f.columns)])
    mypFrame = DataFrame(index=mymarkers.keys(), columns=ids)
    for key, value in mymarkers.items():
        for i in ids:
            mypFrame.loc[key,i] = -10*log(mannwhitneyu(x=f.loc[f['Description'].isin(value),:][i], 
                     y=f[i],alternative='greater').pvalue)
    return(mypFrame)


def add_anno(adata,cnames,mycol,clusters='louvain'):
    """ Adds annotation generated with make_anno to a AnnData object
    Takes as input the AnnData object to which annotation will be appended and the annotation pd
    Generates a pd.Series that can be added as a adata.obs column with annotation at a level
    parameters
    ----------
    adata: AnnData
      AnnData object that is to be annotated
    cnames: panda.DataFrame
      a list of cluster names generated with make_anno
    mycol: string
      column to be added
    cluster: string
      initial cluster used for annotation (column in adata.obs)

    returns
    -------
    pd.Series
        the annotated cells
    """

    newlab=adata.obs[clusters]
    #mycol='celltype0'
    for i in list(set(newlab)):
        newlab=newlab.replace(i, cnames.loc[i,mycol]).copy()
    return(newlab)

def add_anno_torem(adata,cNames, cluster='louvain'):
    """ Adds annotation generated with make_anno to a AnnData object
    Takes as input the AnnData object to which annotation will be appended and the large annotation
    Generates 5 new adata.obs columns with annotation at distinct levels
    parameters
    ----------
    adata: AnnData
      AnnData object that is to be annotated
    cNames: panda.DataFrame
      a list of cluster names generated with make_anno
    cluster: string
      initial cluster used for annotation (column in adata.obs)

    returns
    -------
    AnnData
        the annotated adata object
    """
    adata.obs['clusterID'] = adata.obs[cluster]
    adata.rename_categories('clusterID', cNames)
    
    ### split according to T, B, monocytes, tumor cells
    cellnames = []
    cellgroups = []
    scellgroups = []
    sscellgroups = []
    for i in adata.obs['clusterID']:
        cellnames.append(i.split(".")[0]+i.split(".")[2]+i.split(".")[3])
        cellgroups.append(i.split(".")[2])
        scellgroups.append(i.split(".")[2]+i.split(".")[3])
        sscellgroups.append(i.split(".")[2]+i.split(".")[3]+i.split(".")[4])
    adata.obs['cell_names'] = cellnames
    adata.obs['cell_group'] = cellgroups
    adata.obs['scell_group'] = scellgroups
    adata.obs['sscell_group'] = sscellgroups
    return(adata)


def read_annotconfig(configfile):
    """ Reads the configuration file for signature-based hierarhical annotation. 
    parameters
    ----------
     configfile: 'str'
      path to the configuration file, tsv

    returns
    -------
    sigconfig: panda.DataFrame
      a dataframe with the config details
    levsk: list of str
        String of distinct levels based on configuration file
    """
    ### read the config file
    sigconfig=read_csv(configfile, sep='\t', index_col=0)
    ### Reorder with the specified order. Place better signatures first, only first match will be kept
    sigconfig=sigconfig.sort_values('Order') 
    #### Consider up to 7 levels 
    nochild=list(set(sigconfig.index)-set(sigconfig['Parent']))
    levs=[]
    levs.append(list(sigconfig.loc[sigconfig['Parent']=='None',].index))
    for i in range(0,7):
        levs.append(list(sigconfig.loc[sigconfig['Parent'].isin(levs[i]),].index))

    levsk=[]
    for lev in levs: 
        if len(lev)>0: levsk.append(lev)
    
    return(sigconfig,levsk)


def make_anno(df,sigscores,sigconfig,levsk,lab='celltype'):
    """ Annotate cell types
    Based on a dataframe of -log10pvals, a cutoff and a signature set generate cell annotation
    Hierarchical model of Immune cell annotation.
    It expects a specific set of signatures with a common prefix (signame) and specified suffixes indicating
    the cell type.
    parameters
    ----------
    mypFrame: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    sigscores: dict
      a dictionary with cluster attribution per signature
    sigconfig: panda.DataFrame
      a dataframe with the configuration information
    levsk: list of str
      String of distinct levels based on configuration file
    lab: 'str'
      cell type category base name
      
    returns
    -------
    cnames: panda.DataFrame
      a dataframe with cluster to cell type attribution per distinct levels
    """

    #### First part, get cluster identities
    myclust=list(df.columns)
    annol={}
    for ii in myclust:
        annol[ii]=[]
        for lev in levsk[0]:
            if ii in sigscores[lev]:
                annol[ii]=lev
                break
        if (len(levsk)>2):
            for j in range(len(levsk)-1):
                #jj=j+1
                if len(annol[ii])==0:
                    break
                if len(str.split(annol[ii],"."))< (j+1):
                    break
                sublev=list(sigconfig.loc[sigconfig["Parent"]==str.split(annol[ii],".")[j],:].index)
                sublevk=[]
                for x in sublev:
                    if x in set(list(sigscores.keys())): sublevk.append(x)
                sublev=sublevk.copy()
                for lev in sublev:
                    if ii in sigscores[lev]:
                        annol[ii]=annol[ii]+'.'+lev
                        break
                        
    #### Second part, transform to pandas
    cnames=[]
    for key, val in annol.items():
        if len(val)==0:
            tmp=['Cell']
        else:
            tmp=str.split(val,".")
        for i in range(len(levsk)-len(tmp)):
            tmp.append(tmp[len(tmp)-1])
        cnames.append(tmp)
    cnames=DataFrame(cnames)
    cnames.index=annol.keys()
    cnames.columns=[lab+str(x) for x in list(cnames.columns)]
    
    return(cnames)

def make_anno_torem(mypFrame, myc, signame, f, CD45threshold=0.3,species='human'):
    """ Annotate Immune cells and some other cell types (melanoma, endothelial)
    Based on a dataframe of -log10pvals, a cutoff and a signature set generate cell annotation
    Hierarchical model of Immune cell annotation.
    It expects a specific set of signatures with a common prefix (signame) and specified suffixes indicating
    the cell type.
    parameters
    ----------
    mypFrame: panda.DataFrame
      a dataframe of -log10pvals per signature and cell cluster
    myc: 'numpy.float64'
      threshold used for cluster attribution
    f: panda.DataFrame
      a dataframe with fraction genes expressed per cluster
    signame: 'str'
      signature base name

    returns
    -------
    list of str
        String of cluster labels
    """
    ## over CD45threshold % of cells have CD45
    aCD45P = set(where(gene_fr(f, ['PTPRC']) >= CD45threshold)[0])
    if species == 'mouse':
        aCD45P = set(where(gene_fr(f, ['Ptprc']) >= CD45threshold)[0])
    aBcells = getset(mypFrame, signame+'_Bcells',myc)
    aPlasma = getset(mypFrame, signame+'_Plasma',myc)
    aTcd8 = getset(mypFrame, signame+'_Tcd8',myc)
    aTcd4 = getset(mypFrame, signame+'_Tcd4',myc*1/2)
    aTcgd = getset(mypFrame, signame+'_Tcgd',myc*2/3) 
    aTreg = getset(mypFrame, signame+'_Treg',myc)
    aTNK = getset(mypFrame, signame+'_TNK',myc*2)
    aTilCM = getset(mypFrame, signame+'_TilCM',myc*2)
    aT4CM = getset(mypFrame, signame+'_T4CM',myc*3/2)
    aTcytox = getset(mypFrame, signame+'_Tcytox',myc*2)
    aTEM = getset(mypFrame, signame+'_TEM',myc*2)
    aTtexh = getset(mypFrame, signame+'_Ttexh',myc)
    aTpexh = getset(mypFrame, signame+'_Tpexh',myc)
    aNKT = getset(mypFrame, signame+'_NKT',myc*3/2)
    aTcells = getset(mypFrame, signame+'_Tcells',myc*4/3)
    aNKcells = getset(mypFrame, signame+'_NKcells',myc*2)
    aNKnai = getset(mypFrame, signame+'_NKnai',myc)
    aNKcyt = getset(mypFrame, signame+'_NKcyt',myc*2)
    aMyelo = getset(mypFrame, signame+'_Myelo',myc)
    aNai = getset(mypFrame, signame+'_Naive',myc)
    aAct = getset(mypFrame, signame+'_Activation',myc)
    aMem = getset(mypFrame, signame+'_Memory',myc)
    aEff = getset(mypFrame, signame+'_Eff',myc)
    aNonEff = getset(mypFrame, signame+'_NonEff',myc)
    aCytotox = getset(mypFrame, signame+'_Cytotox',myc*2)
    aDCs = getset(mypFrame, signame+'_aDCs', myc)
    moDC = getset(mypFrame, signame+'_moDC', myc*3/2) 
    pDCs = getset(mypFrame, signame+'_pDCs', myc*2)   
    acDC1 = getset(mypFrame, signame+'_cDC1', myc)       
    acDC2 = getset(mypFrame, signame+'_cDC2', myc*2)       
    aMono = getset(mypFrame, signame+'_Monocytes', myc)
    aTAM = getset(mypFrame, signame+'_TAM', myc*3/2)   
    aTMO = getset(mypFrame, signame+'_TMO', myc)
    aMo14 = getset(mypFrame, signame+'_Mo14', myc)
    aMo16 = getset(mypFrame, signame+'_Mo16', myc)
    aTAMCx = getset(mypFrame, signame+'_TAMCx', myc)
    aTMid = getset(mypFrame, signame+'_TMid', myc)
    aMoMa = getset(mypFrame, signame+'_MoMa', myc*2)
    aMyeloSubtype = getset(mypFrame, signame+'_MyeloSubtype', myc)
    aGranu = getset(mypFrame, signame+'_Granulo', myc)
    aNeutro = getset(mypFrame, signame+'_Neutrophil', myc)
    aMacrophage = getset(mypFrame, signame+'_Macrophage', myc)
    aCC = getset(mypFrame, signame+'_Cellcycle', myc*4/3)
    aCh = getset(mypFrame, signame+'_Checkpoint', myc)
    aIfng = getset(mypFrame, signame+'_Ifng', myc)
    aEndo = getset(mypFrame, signame+'_Endo', myc)
    aCafs = getset(mypFrame, signame+'_Cafs', myc)
    aMelMelan = getset(mypFrame, signame+'_MelMelan', myc)
    aMega = getset(mypFrame, signame+'_Megakaryocytes', myc)
   #### Part 2 Combine
    CD45P = set(aCD45P).union(set(aBcells).union(aNKcells).union(aTNK).union(aTcells).union(aTreg).union(aMyelo).union(pDCs))
    Tc = set(CD45P).intersection(set(aTcells).union(aTreg)-set(aBcells)-set(aMyelo))
    Nk = set(CD45P).intersection(set(aNKcells).union(set(aCytotox))-set(aTcells))
    Pl = set(aPlasma).intersection(set(CD45P))-Tc-Nk
    Bc = set(aBcells).union(Pl)
    Tr = set(Tc).intersection(aTreg)
    Dc = set(aDCs).union(acDC1).union(acDC2).union(moDC).union(pDCs)-Tc-Nk-Bc-Pl
    My = set(CD45P).intersection(set(aMoMa).union(set(aTAM).union(set(aTMO)).union(set(aTMid)))).union(Dc)-Tc-set(aBcells)-Nk
    naiTc = set(Tc).intersection(set(aNai))
    ccTc = set(Tc).intersection(set(aCC))
    chTc = set(Tc).intersection(set(aCh))
    effTc = set(Tc).intersection(set(aCytotox)-Tr)
    actTc = set(Tc).intersection(set(aAct)-Tr)
    memTc = set(Tc).intersection(set(aMem))-actTc-Tr-naiTc
    exTc = set(Tc).intersection(chTc.union(aNonEff))-naiTc
    texTc = set(Tc).intersection(aTtexh)-naiTc
    pexTc = set(Tc).intersection(aTpexh.intersection(chTc))-naiTc
    cytoTc = set(Tc).intersection(set(aTcytox)-Tr)-naiTc
    NKTc = set(Tc).intersection(set(aNKT)-Tr)
    TilCM = set(Tc).intersection(set(aTilCM)-Tr) #IL7R cluster
    T4CM = set(Tc).intersection(set(aT4CM)-Tr)-TilCM-cytoTc
    TEM = set(Tc).intersection(set(aTEM).intersection(aTcytox)-Tr).union(aEff)-TilCM-T4CM-cytoTc
    TcCD8 = set(aTcd8).intersection(Tc)
    Tr = Tr-TcCD8
    TcCD4 = set(aTcd4).intersection(Tc)-Tr-set(TcCD8)
    Tgd = set(Tc).intersection(aTcgd)
    Nknai = aNKnai.intersection(Nk)
    Nkcyt = aNKcyt.intersection(Nk)-Nknai
    Mo14 = aMo14.intersection(My)-acDC1-acDC2-aDCs-Tc-Bc-Nk
    Mo16 = aMo16.intersection(My)-acDC1-acDC2-aDCs-Tc-Bc-Nk
    Mo14 = Mo14-Mo16
    pDCs = pDCs-Bc
    #### Part 3 annotate
    ids = list(set(f.iloc[:, 2:len(f.columns)]))
    nclust = len(ids)
    ###### Summarize the annotations:
    cNames = []
    CD45P = [str(i) for i in CD45P]
    for i in range(0,nclust):
        ii = str(i)
        newname = "C"+ii+['.nCD45', '.CD45'][(ii in CD45P)*1]
        if (ii not in CD45P):
            lev2 = ['no', '.Endo'][(ii in aEndo)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.Mel'][(ii in aMelMelan)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.Caf'][(ii in aCafs)*1]
            lev3='.'
            lev4='.'
            
        if (ii in CD45P):
            lev2 = ['no', '.Tc'][(ii in Tc)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.NK'][(ii in Nk)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.Bc'][(ii in Bc)*1]
            if (lev2 == 'no'):
                lev2 = ['no', '.My'][(ii in My)*1]
            if (lev2 == 'no'):
                lev2='.'
                lev3='.'
                lev4='.'
    
        if (lev2 == '.Tc'):
            lev3 = ['no', '.8'][(ii in TcCD8)*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Tr'][(ii in Tr)*1]
            if (lev3 == 'no'):
                lev3 = ['.4', '.4'][(ii in TcCD4)*1]
            lev4 = ['no', '.Tr'][(ii in Tr)*1]          
            if (lev4 == 'no'):
                lev4 = ['no', '.texh'][(ii in texTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.NKT'][(ii in NKTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Tgd'][(ii in Tgd)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Cytox'][(ii in cytoTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.CMil7'][(ii in TilCM)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.CM4'][(ii in T4CM)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.CC'][(ii in ccTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.texh'][(ii in pexTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.EM'][(ii in TEM)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Nai'][(ii in naiTc)*1]  
            if (lev4 == 'no'):
                lev4 = ['no', '.Exh'][(ii in chTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Eff'][(ii in effTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Act'][(ii in actTc)*1]
            if (lev4 == 'no'):
                lev4 = ['no', '.Mem'][(ii in memTc)*1]
            if (lev4 == 'no'):
                lev4 = '.Mid'
        if (lev2 == '.Bc'):
            lev3 = ['.Bc', '.Pl'][(ii in Pl)*1]
            lev4 = '.'
        if (lev2 == '.NK'):
            lev3 = ['no', '.Nai'][(ii in Nknai)*1]
            if (lev3 == 'no'):
                lev3 = '.Cyt'
            lev4 = '.'               
        if (lev2 == '.My'):
            lev3 = ['no', '.pDC'][(ii in pDCs.intersection(My))*1]   
            if (lev3 == 'no'):
                lev3 = ['no', '.aDCs'][(ii in aDCs.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.cDC1'][(ii in acDC1.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.cDC2'][(ii in acDC2.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.aTc'][(ii in aTcells.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Mo14'][(ii in Mo14.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Mo16'][(ii in Mo16.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.TAMCx'][(ii in aTAMCx.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.TMid'][(ii in aTMid.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.aTAM'][(ii in aTAM.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.moDC'][(ii in moDC.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.TMO'][(ii in aTMO.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.MoMa'][(ii in aMyeloSubtype.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Gra'][(ii in aGranu.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Mky'][(ii in aMega.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Neu'][(ii in aNeutro.intersection(My))*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.Tc'][(ii in aTcells)*1]
            if (lev3 == 'no'):
                lev3 = ['no', '.My'][(ii in aMyelo)*1]
            if (lev3 == 'no'):
                lev3 = '.'
            lev4 = '.'
        newname = newname+lev2+lev3+lev4+['.nNa', '.Na'][(ii in aNai)*1]
        newname += ['.nCC', '.CC'][(ii in aCC)*1]+['.nCy', '.Cy'][(ii in aCytotox)*1]
        newname += ['.nCh', '.Ch'][(ii in aCh)*1]+['.nAc', '.Ac'][(ii in aAct)*1]
        cNames.append(newname)
    return(cNames)

def match_cluster(adata,obsquery,obsqueryval,obsref='leiden',cutoff=0.5):
    """ Matches categories from adata.obs to each other. For a query category specified in obsquery
    and a value specified in obsqueryval, checks which clusters (or other adata.obs categories, obsref)
    contain >50% (or distinct cutoff, cutoff) of cells of the specified kind.
    
    parameters
    ----------
    adata: AnnData
      AnnData object 
    obsquery: 'str'
      adata.obs category name used for querying
    obsqueryval: 'str'
      adata.obs category name value, present in obsquery
    obsref: 'str'
      adata.obs category name to be returned
    cutoff: 'numpy.float64'
      fraction of positive cells returned
    returns
    -------
    list
        a list of the cluster IDs that match the query label
    """
    
    myc=[]
    myleiden=list(set(adata[adata.obs[obsquery].isin([obsqueryval])].obs[obsref]))
    for i in myleiden:
        mysub=adata[adata.obs[obsref].isin([i])].copy()
        myc.append(len(mysub[mysub.obs[obsquery]==obsqueryval])/len(mysub))
    return(list(array(myleiden)[array(myc)>cutoff]))

