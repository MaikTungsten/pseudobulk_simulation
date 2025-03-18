import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from tqdm import tqdm
import random
from numpy.random import choice
import warnings
from scipy.stats import rankdata


# Simulate data

def generate_simulated_data(scData, outname=None,
                            d_prior=None,
                            n=1000, samplenum=5000, props = None,
                            random_state=None, sparse=True, sparse_prob=0.2,
                            rare=False, rare_percentage=0.3):
    
    # scData is a pd.DataFrame (a cells x genes matrix)
    
    genes = list(scData.iloc[:,0:scData.shape[1]-1].columns) # save gene names for later

    # Get number and groups of cell types present
      
    num_celltype = len(scData['CellType'].value_counts())
    celltype_groups = scData.groupby('CellType').groups # dictionary of cell types
    scData.drop(columns='CellType', inplace=True)
    
    if not isinstance(props, pd.DataFrame):
        # generate random cell type proportions

        if random_state is not None and isinstance(random_state, int):
            print('Random state specified. This will improve the reproducibility.')

        if d_prior is None:
            print('Generating cell fractions using Dirichlet distribution without prior info (actually random).')
            if isinstance(random_state, int):
                np.random.seed(random_state)
            prop = np.random.dirichlet(np.ones(num_celltype), samplenum) # randomly alter proportion of 1 per cell type via dirichlet distribution for samplenum, e.g. 5000 times
            print('Random cell fractions are generated.')
        elif d_prior is not None:
            print('Using prior info to generate cell fractions in Dirichlet distribution')
            assert len(d_prior) == num_celltype, 'dirichlet prior is a vector, its length should equals ' \
                                                'to the number of cell types'
            if isinstance(random_state, int):
                np.random.seed(random_state)
            prop = np.random.dirichlet(d_prior, samplenum) # d_prior would set an initial cell type distribution to randomly and slightly alter for samplenum
            print('Dirichlet cell fractions are generated.')

        prop = prop / np.sum(prop, axis=1).reshape(-1, 1) # scale to have percentage proportions for each cell type summing up to 100 %; each row is one simulated sample and its proportions
        
        
        # sparse cell fractions
        if sparse:
            print("You set sparse as True, some cell's fraction will be zero, the probability is", sparse_prob)
            ## Only partial simulated data is composed of sparse celltype distribution
            for i in range(int(prop.shape[0] * sparse_prob)):
                indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * sparse_prob)) # chose a random cell type index from 0 to no. of cell types to drop out per sample to be simulated
                prop[i, indices] = 0

            prop = prop / np.sum(prop, axis=1).reshape(-1, 1) # resize proportions

        if rare:
            print(
                'Selected rare, thus some cell type fractions are very small (<3%), '
                'this celltype is randomly chosen by percentage you set before.')
            ## choose celltype to be rare
            np.random.seed(0)
            indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * rare_percentage))
            prop = prop / np.sum(prop, axis=1).reshape(-1, 1)

            for i in range(int(0.5 * prop.shape[0]) + int(int(rare_percentage * 0.5 * prop.shape[0]))):
                prop[i, indices] = np.random.uniform(0, 0.03, len(indices)) # rare between 0 % and 3 %
                buf = prop[i, indices].copy()
                prop[i, indices] = 0
                prop[i] = (1 - np.sum(buf)) * prop[i] / np.sum(prop[i])
                prop[i, indices] = buf
    
    elif isinstance(props, pd.DataFrame):
        prop = np.array(props)

    # precise number for each celltype
    cell_num = np.floor(n * prop) # from proportions get the number of cells per cell type to be sampled given n (e.g. n=1000 cells)

    # precise proportion based on cell_num
    prop = cell_num / np.sum(cell_num, axis=1).reshape(-1, 1) # then, obtain accurate proportions that fit the given number n

    # start sampling
    sample = np.zeros((prop.shape[0], scData.shape[1])) # zeros for samplenum x no. embeddings
    allcellname = celltype_groups.keys() # cell type names
    print('Sampling cells to compose pseudo-bulk data...')
    
    for i, sample_prop in tqdm(enumerate(cell_num)): # tqdm shows nice progress meter, enumerate lists line number and number of cells for each cell type to be sampled
        for j, cellname in enumerate(allcellname):
            select_index = choice(celltype_groups[cellname], size=int(sample_prop[j]), replace=True) # choose cell type names for a given number of cells to be sampled for the type
            sample[i] += (scData.loc[select_index,:]).sum(axis=0) # linear combination of cells to obtain aggregared read count per gene across sampled cells and cell types

    if not isinstance(props, pd.DataFrame):
        sampleDF = pd.DataFrame(sample, index=['Sample_'+str(i) for i in range(1,prop.shape[0]+1)], columns = genes[0:sample.shape[1]])
        prop = pd.DataFrame(prop, index = ['Sample_'+str(i) for i in range(1,prop.shape[0]+1)], columns=celltype_groups.keys())

    elif isinstance(props, pd.DataFrame):
        sampleDF = pd.DataFrame(sample, index = props.index.tolist(), columns = genes[0:sample.shape[1]])
        prop = pd.DataFrame(prop, index = props.index.tolist(), columns = props.columns.tolist())
    
    simudata = anndata.AnnData(X=sampleDF,
                               obs=prop,
                               #var=genes[0:sample.shape[1]]
                               )

    print('Sampling is done.')
    #if outname is not None:
    #    simudata.write_h5ad(outname + '.h5ad')
    
    return simudata # anndata object containing ground-truth proportions in observations and having geneIDs as numbers for hidden layers



# Select mRNA genes

def main_gene_selection(X_df, gene_list):
    """
    Describe:
        rebuild the input adata to select target genes that encode proteins 
    Parameters:
        adata->`~anndata.AnnData` object: adata with var index_name by gene symbol
        gene_list->list: wanted target gene 
    Returns:
        adata_new->`~anndata.AnnData` object
        to_fill_columns->list: zero padding gene
    """
    X_df = X_df.fillna(0)
    to_fill_columns = list(set(gene_list) - set(X_df.columns))

    padding_df = pd.DataFrame(np.zeros((X_df.shape[0], len(to_fill_columns))), 
                              columns=to_fill_columns, 
                              index=X_df.index)
    
    X_df = pd.DataFrame(np.concatenate([df.values for df in [X_df, padding_df]], axis=1), 
                        index=X_df.index, 
                        columns=list(X_df.columns) + list(padding_df.columns))
    
    X_df = X_df[gene_list]
    
    var = pd.DataFrame(index=X_df.columns)
    var['mask'] = [1 if i in to_fill_columns else 0 for i in list(var.index)]
    
    return X_df, to_fill_columns, var

def normRank(pseudobulk):
    df = np.array(pseudobulk.to_df())

    ranked_array = rankdata(df, axis = 1, method = 'min')
    ranked_array = ranked_array/ranked_array.shape[1]

    ranked_df = pd.DataFrame(ranked_array)
    ranked_df.index = pseudobulk.obs_names.to_list()
    ranked_df.columns = pseudobulk.var_names.to_list()

    pseudobulk.layers['ranked'] = ranked_df
    return pseudobulk
