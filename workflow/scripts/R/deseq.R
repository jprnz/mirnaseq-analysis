log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Input
countsfile <- snakemake@input[['counts']]
samplefile <- snakemake@input[['samples']]
configfile <- snakemake@input[['config']]
analysis_name <- snakemake@wildcards[['analysis']]
dataset_name <- snakemake@wildcards[['dataset']]
threads <- snakemake@threads

# Output
output_xls <- snakemake@output[['xls']]
output_cnt <- snakemake@output[['cnt']]
output_norm <- snakemake@output[['norm']]

library("DESeq2")
library("RColorBrewer")
library("WriteXLS")
library("data.table")
library("pheatmap")
library("yaml")
library("ggplot2")

# Set seed, scipen and load functions
source(file.path(snakemake@scriptdir, "functions.R"))

# Create places for supplementary outputs
prefix <- dirname(output_xls)
rnkpath <- file.path(prefix, "rnk-files")
tabpath <- file.path(prefix, "tables")
dir.create(rnkpath, recursive=TRUE, showWarnings=FALSE)
dir.create(tabpath, recursive=TRUE, showWarnings=FALSE)

# Read config
config <- read_config(configfile, analysis_name)

# Read coldata
coldat <- read_samplefile(samplefile)

# Read counts for samples in coldat
count <- read_counts(countsfile, coldat)
count <- count[apply(count, 1, max) >= 10, ]

# Get model
model <- as.formula(config$design)

# Read genes file
#genes <- fread(genesfile)

# Make DESeq data object and run analysis
dds <- DESeqDataSetFromMatrix(countData=count, colData=coldat, design=model)
dds <- do.call(DESeq, c(list(dds), config$deseq))

# Get comparisons
cat("ResultsNames:\n ", paste0(resultsNames(dds), collapse="\n  "), "\n")
comparisons <- get_comparisons(config$results, dds)

# Get results
res = get_results(dds, comparisons, config$results)

# Make output tables and rnk files
tab <- list()
grp <- get_group_counts(dds, config)
for (cmp in names(res)) {
    tab[[cmp]] <- format_tab(res[[cmp]], dds, grp)
}

# Write results
for (cmp in names(res)) {
   tab_filename <- paste0(tabpath, "/", cmp, ".csv")
   fwrite(tab[[cmp]], file=tab_filename)

   plt_filename <- paste0(prefix, "/Heatmap-", cmp, ".pdf")
   plot_heatmap(res, cmp, dds, config, comparisons, plt_filename)

   plt_filename <- paste0(prefix, "/Volcano-", cmp, ".pdf")
   plot_volcano(res[[cmp]], plt_filename)
}

# Plot disp estimates
plt_filename <- paste0(prefix, "/Dispersions.pdf")
plot_dispersions(dds, plt_filename)

# Write PCA plot
config$pca$filename = paste0(prefix, "/PCA.pdf")
do.call(plot_pca, c(list(counts(dds, normalized=TRUE), coldat), config$pca))

# Write Excel
cat("Writting XLS...\n")
write_xls(tab, comparisons, filename=output_xls)

# Write counts
cat("Writting Counts...\n")
write_counts(counts(dds), output_cnt)
write_counts(counts(dds, normalized=TRUE), output_norm)

# Write workspace
cat("Writting data image...\n")
save.image(file=paste0(prefix, "/analysis.Rdata"))

