# This documents two command line applications for pseudobulk simulations

For pseudobulk simulations, the Scaden-based pseudobulk simulator slightly adapted and the SimBu-based simulator can be used. For each case, a little command-line application is available.

## Scaden-based pseudobulk simulation

The Scaden-based simulator is used to simulate a variety of pseudobulks with diverse cell type proportions. With the sparse and rare parameters, the simulation can be fine-tuned to randomly set cell type proportions to 0 in a specific pseudobulk. On this basis, cells will be sampled from the scRNA-seq data and aggregated, i.e. summed, and scaled to counts per million (CPM), normalized by rank or kept in raw counts scale. You can specify the number of cells to use to build each pseudobulk and the number of pseudobulks to simulate. Importantly, the scRNA-seq data should be in h5ad format and contain a column named "cell_type" in the observations. Based on this column, cell types are identified, proportions simulated and pseudobulks constructed. The returned pseudobulks will be in csv format.

```
usage: simulation.py [-h] [--sc_path SC_PATH] [--sc_layer SC_LAYER] [--samplenum SAMPLENUM] [--props PROPS] [--ncells NCELLS] [--rare RARE] [--norm {rank,CPM,raw}] [--filter_genes {all, mRNA, 3k}] [--sparse SPARSE] [--outname OUTNAME]

Basic tool for simulating pseudobulk data from single-cell data.

options:
  -h, --help            show this help message and exit
  --sc_path SC_PATH     Input single-cell data path.
  --sc_layer SC_LAYER   Input layer of single-cell data to use for normalization.
  --samplenum SAMPLENUM
                        Number of pseudobulk samples to be simulated.
  --props PROPS         If desired, specify path to a csv file that provides proportions for the desired amount of samples. Must contain all cell types given in sc data.
  --ncells NCELLS       Number of pseudobulk samples to be simulated.
  --rare RARE           Probability for rare cell types.
  --norm {rank,CPM,raw}
                        Normalization strategy for pseudobulk after aggregation of single cells.
  --filter_genes {all, mRNA, 3k}
                        Selection of set of genes from pseudobulks: mRNA, all genes as in single-cels, 3k most variable genes.
  --sparse SPARSE       Probability for sparse cell types (e.g. cell type is not present).
  --outname OUTNAME     Ideally, specifiy a date or tissue ID.
```

## SimBu-based pseudobulk simulation

While SimBu is very versatile itself regarding the options for cell type proportions (e.g. equally, randomly distributed or specified). However, here the user needs to provide cell type proportions for a desired amount of samples to be simulated. If randomly simulated proportions are desired, the Scaden-based simulator output of cell type proportions may be a suitable choice. In addition, the output folder and the input scRNA-seq data (in h5ad format) need to be specified. Make sure that cell types are properly annotated in the obs of your h5ad object and that they match the column names of your cell type proportions. The filtering option is currently not adapted to a count cut-off of 3000, as this may remove cell types from your input data that are contained in your specified proportions. Filtering only removes genes for which expression of 0 is recorded for all cells in the scRNA-seq data.

Run the simulation in the given miniconda environment "pseudobulksSimBu" using this command:
```
Rscript Simulation_with_SimBu.R [-h] [-p PROP_PATH] [-c SC_PATH] [-o OUT_PATH]
                               [-f FILTERING]

options:
  -h, --help            show this help message and exit
  -p PROP_PATH, --prop_path PROP_PATH
                        Path to cell type proportions for simulation of
                        pseudobulk data.
  -c SC_PATH, --sc_path SC_PATH
                        Path to single-cell sequencing data for deconvolution.
                        Data format should be cells x genes, which are used as
                        row and colnames.
  -o OUT_PATH, --out_path OUT_PATH
                        Output path for analysis results.
  -f FILTERING, --filtering FILTERING
                        If T, filtering will remove cells less abundant than
                        3000 cells in the scRNA-seq data and remove zero-
                        expressed genes.
```
The generated pseudobulks will be scaled to counts per million (CPM) using two different inputs for simulation: raw read counts (Output file: Pseudobulk_counts_scaled.csv) and CPM-normalized counts (Output file: Pseudobulk_counts_cpm.csv).

## Performance of simulators

Without specification of mRNA biases - only possible with SimBu so far - pseudobulks simulated with the Scaden- and SimBu-based techniques are very similar, when using the same cell type proportions. You can find a comparison in the Jupyter notebook.