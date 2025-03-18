### Using SimBu for pseudobulk simulation ###

# Import packages required for computation

print('Loading packages...')
suppressWarnings(suppressMessages(library(argparse)))
suppressWarnings(suppressMessages(library(SingleCellExperiment)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(SeuratObject)))
suppressWarnings(suppressMessages(library(SimBu)))
suppressWarnings(suppressMessages(library(zellkonverter)))


# import arguments from command line

print('Importing arguments...')
parser <- ArgumentParser()

parser$add_argument("-p", "--prop_path", type = "character", default = "./Proportions.csv", 
  help = "Path to cell type proportions for simulation of pseudobulk data.")

parser$add_argument("-c", "--sc_path", default = "Path/to/scData.h5ad", type = "character",
  help = "Path to single-cell sequencing data for deconvolution. Data format should be cells x genes, which are used as row and colnames.")

parser$add_argument("-o", "--out_path", default = "./output_SimBu/", type = "character",
  help = "Output path for analysis results.")

parser$add_argument("-f", "--filtering", default = "F", type = "character",
  help = "If T, filtering will remove cells less abundant than 3000 cells in the scRNA-seq data and remove zero-expressed genes.")

args <- parser$parse_args()

print('Reading h5ad file...')
# Read H5AD file
scData <- readH5AD(args$sc_path)

# Get counts as matrix
scMat <- as.matrix(assay(scData, "X"))

print('Creating seurat object...')
# Create a seurat object
seurat_obj <- CreateSeuratObject(counts = scMat)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "RC", scale.factor = 1e6) # add layer with CPM-normalized data

print('Extracting counts and cpm...')
# Get raw count and cpm normalized data
counts <- as.matrix(seurat_obj@assays$RNA$counts)
cpm <- as.matrix(seurat_obj@assays$RNA$data)

# Get cell type annotation
annotation <- as.data.frame(colData(scData))
annotation$ID <- rownames(annotation)

# Construct the SimBu dataset

print('Constructing SimBu dataset...')
if (args$filtering == "T"){
  ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  tpm_matrix = cpm,
  name = "input_dataset",
  #type_abundance_cutoff = 3000, # remove cell types less abundant than 3000 cells
  filter_genes = TRUE # remove genes with 0 expression in all cells
  )
} else {
  ds <- SimBu::dataset(
  annotation = annotation,
  count_matrix = counts,
  tpm_matrix = cpm,
  name = "input_dataset",
  filter_genes = FALSE # no pre-filtering
  )
}

# Load cell type proportions to specifiy simulation

proportions <- as.data.frame(read.csv(args$prop_path, header = TRUE))
rownames(proportions) <- proportions$X
proportions <- proportions[,2:ncol(proportions)]
colnames(proportions) <- gsub("[.]", " ", colnames(proportions))

# Perform simulation with SimBu

print('Simulating bulk data with given proportions...')
simulation <- SimBu::simulate_bulk(
  data = ds,
  scenario = "custom",
  custom_scenario_data = proportions,
  scaling_factor = "NONE",
  ncells = 1000,
  nsamples = nrow(proportions),
  BPPARAM = BiocParallel::MulticoreParam(workers = 4), # this will use 4 threads to run the simulation
  run_parallel = TRUE,
  total_read_counts = 1e6
)

# Extract simulated data and export

print('Extracting scaled pseudobulk count data...')
# Count input
df <- t(as.data.frame(SummarizedExperiment::assays(simulation$bulk)[["bulk_counts"]]))
bulkName <- paste0(args$out_path, '/Pseudobulk_counts_scaled.csv')
write.csv(df, bulkName)

print('Extracting cpm-normalized pseudobulk data generated from cpm-normalized cells...')
# CPM input
df <- t(as.data.frame(SummarizedExperiment::assays(simulation$bulk)[["bulk_tpm"]]))
bulkName <- paste0(args$out_path, '/Pseudobulk_counts_cpm.csv')
write.csv(df, bulkName)
