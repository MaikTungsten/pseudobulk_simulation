import argparse
import anndata
import numpy as np
import scanpy as sc
import pandas as pd
from tqdm import tqdm
import random
from numpy.random import choice
import warnings
from tools import generate_simulated_data, main_gene_selection, normRank



# Setup argparser
parser = argparse.ArgumentParser(description='Basic tool for simulating pseudobulk data from single-cell data.')

parser.add_argument('--sc_path', type=str, default='./', help='Input single-cell data path.')
parser.add_argument('--sc_layer', type=str, default='unspecified', help='Input layer of single-cell data to use for normalization.')
parser.add_argument('--samplenum', type=int, default=5000, help='Number of pseudobulk samples to be simulated.')
parser.add_argument('--props', type=str, default=None, help='If desired, specify path to a csv file that provides proportions for the desired amount of samples. Must contain all cell types given in sc data.')
parser.add_argument('--ncells', type=int, default=1000, help='Number of pseudobulk samples to be simulated.')
parser.add_argument('--rare', type=float, default=0.3, help='Probability for rare cell types.')
parser.add_argument('--norm', type=str, default='CPM', choices=['rank', 'CPM', 'raw'], help='Normalization strategy for pseudobulk after aggregation of single cells.')
parser.add_argument('--filter_genes', type=str, default='mRNA', choices=['all', 'mRNA', '3k'], help='Selection of set of genes from pseudobulks: mRNA, all genes as in single-cels, 3k most variable genes.')
parser.add_argument('--sparse', type=float, default=0.2, help='Probability for sparse cell types (e.g. cell type is not present).')
parser.add_argument('--outname', type=str, default='NoDate', help='Ideally, specifiy a date or tissue ID.')

args = parser.parse_args()

# Extract arguments
ncells = args.ncells
samplenum = args.samplenum
sc_path = args.sc_path
outputname = args.outname
sparse = args.sparse
rare = args.rare
propPath = args.props
sc_layer = args.sc_layer
norm = args.norm
filter_genes = args.filter_genes

# Step 1: Load single-cell data, select relevant pre-normalited layer; scData should already be QC'ed beforehand; load proportions if supplied

# Load single-cell data
inData = sc.read_h5ad(sc_path)

if sc_layer not in inData.layers.keys():
    if sc_layer == "unspecified":
        sc_layer = 'X'
    else:
        print('Specified sc_layer ' + sc_layer + ' not in scRNA-seq data object. Using default layer.')
        sc_layer = 'X'
        

# Extract relevant layer
if sc_layer == 'X':
    print('Using default layer for pseudobulk simulation.')
    scData = inData.to_df() # get X to dataframe
    scData['CellType'] = inData.obs['cell_type'] # annotate a CellType column from cell_type column in obs
elif sc_layer != 'X':
    print('Using layer ' + sc_layer + ' from scRNA-seq data for pseudobulk simulation.')
    scData = inData.to_df(layer = sc_layer) # get X to dataframe
    scData['CellType'] = inData.obs['cell_type'] # annotate a CellType column from cell_type column in obs

# Load proportions, if supplied

if propPath == None:
    props = None
elif propPath.endswith('.csv'):
    props = pd.read_csv(propPath, index_col=0)
    if props.shape[0] != samplenum:
        print('Number of samples in proportions is not matching specified sample number. Sample number is now adjusted to ' + str(props.shape[0]) +'.')
        samplenum = props.shape[0]
    elif props.shape[1] != len(scData['CellType'].value_counts()):
        raise ValueError('Number of cell types in proportions is not matching number of cell types in scRNA-seq data.')
    else:
        print('Props match specified parameters.')
else:
    raise ValueError('Proportions not in csv format. Please supply as csv file of samples x cell types.')


# Step 2: Create simulated data

pseudobulks = generate_simulated_data(scData,
                                      n = ncells,
                                      samplenum = samplenum,
                                      props=props,
                                      sparse=True,
                                      sparse_prob=sparse,
                                      rare=True,
                                      rare_percentage=rare)


# Step 3: Normalize simulated data and filter data in different scenarios

# CPM normalization option
if norm == 'CPM':
    print('Scaling pseudobulks to CPM.')
    sc.pp.normalize_total(pseudobulks, target_sum=1e6)

    ##### All genes #####
    if filter_genes == "all":
        pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions
        
    ##### Only mRNA genes #####
    elif filter_genes == "mRNA":
        pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

        # Import gene list for filtering
        gene_list_df = pd.read_csv('mRNA_annotation.tsv', header=0, delimiter='\t')
        gene_list = list(gene_list_df['gene_name'])

        # Select (and add) genes as necessary
        pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)

    ##### 3000 highly variable genes; not scaled #####
    elif filter_genes == "3k":
        highlyVariableDF = sc.pp.log1p(pseudobulks, copy = True)
        x = sc.pp.highly_variable_genes(highlyVariableDF, n_top_genes=3000, inplace=False)
        y = x.loc[x['highly_variable']==True].index.to_list()
        pseudobulkDF = pd.DataFrame(pseudobulks[:,y].X, index=pseudobulks[:,y].obs_names, columns=pseudobulks[:,y].var_names)
        proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

# Rank normalization option
elif norm == 'rank':
    print('Ranking genes in pseudobulks.')
    pseudobulks = normRank(pseudobulks)

    pseudobulkDF = pd.DataFrame(pseudobulks.layers['ranked'], index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
    proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

    ##### Only mRNA genes #####
    if filter_genes == "mRNA":
        # Import gene list for filtering
        gene_list_df = pd.read_csv('mRNA_annotation.tsv', header=0, delimiter='\t')
        gene_list = list(gene_list_df['gene_name'])

        # Select (and add) genes as necessary
        pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)

    ##### 3000 highly variable genes; not scaled #####
    elif filter_genes == "3k":
        highlyVariableDF = sc.pp.log1p(pseudobulks, copy = True)
        x = sc.pp.highly_variable_genes(highlyVariableDF, n_top_genes=3000, inplace=False)
        y = x.loc[x['highly_variable']==True].index.to_list()
        pseudobulkDF = pd.DataFrame(pseudobulks[:,y].layers['ranked'], index=pseudobulks[:,y].obs_names, columns=pseudobulks[:,y].var_names)
            
# Raw counts option - no normalization
elif norm == 'raw':
    print('Returning raw summed counts in pseudobulks.')
    pseudobulkDF = pd.DataFrame(pseudobulks.X, index=pseudobulks.obs_names, columns=pseudobulks.var_names) 
    proportionsDF = pd.DataFrame(pseudobulks.obs) # does not need normalization --> proportions

    ##### Only mRNA genes #####
    if filter_genes == "mRNA":
        # Import gene list for filtering
        gene_list_df = pd.read_csv('mRNA_annotation.tsv', header=0, delimiter='\t')
        gene_list = list(gene_list_df['gene_name'])

        # Select (and add) genes as necessary
        pseudobulkDF, _, _ = main_gene_selection(pseudobulkDF, gene_list)

    ##### 3000 highly variable genes; not scaled #####
    elif filter_genes == "3k":
        highlyVariableDF = sc.pp.log1p(pseudobulks, copy = True)
        x = sc.pp.highly_variable_genes(highlyVariableDF, n_top_genes=3000, inplace=False)
        y = x.loc[x['highly_variable']==True].index.to_list()
        pseudobulkDF = pd.DataFrame(pseudobulks[:,y].X, index=pseudobulks[:,y].obs_names, columns=pseudobulks[:,y].var_names)


#Step 4: Export as csv
print('Writing pseudobulks to output.')

proportionsDF.to_csv(outputname+'_pseudobulk_proprotions.csv')
pseudobulkDF.to_csv(outputname+'_pseudobulks.csv')