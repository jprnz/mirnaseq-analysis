# miRNA-Seq Workflow
Mapping and Differential expression of miRNA-Seq data

This workflow is designed to perform simple case-control analysis of one or more variables using DESeq.
The workflow leverages ties together the following pieces of software:
* [fastp](https://github.com/OpenGene/fastp)
* [umi-tools](https://umi-tools.readthedocs.io/en/latest)
* [umicollapse](https://github.com/Daniel-Liu-c0deb0t/UMICollapse)
* [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
* [MultiQC](https://multiqc.info/)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

All software dependencies are automatically resolved using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/index.html).
Local execution and running via [SLURM](https://slurm.schedmd.com/) are supported.  
When running from an HPC all work should be done on an interactive session.

# Setup workflow
## Code
Clone this repository into your project directory
```sh
git clone <url>
```

## Conda
Install conda and activate the snakemake environment
```sh
./setup
source conda/bin/activate snakemake
```

# Running
See Snakemake documentation for more [command line options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options).
Default command line options are set using `--workflow-profiles` and can be found [here](workflow/profiles/config.yaml).
You can override these pre-sets by using their command-line equivalents or clear these using `--workflow-profiles none`

### Cluster usage
When working on a cluster it is generally best to be running snakemake in an interactive session
```sh 
srun --mem 8G --pty bash
```

### Run simple test (optional)
To run the example data using SLURM, from the `.test/` directory, run
```sh
cd .test
snakemake -s ../workflow/Snakefile --slurm
```

## Running
Setup genome and ensure the settings in `config/config.yaml` are correct, then prepare `fastqs.tsv`, `samples.tsv`, `analysis.yaml`.
See documentation [here](config/README.md) for more information.

The entire workflow can be run without additional options (in this example using SLURM):
```sh
snakemake --executor slurm -A <account>
```

## Pushing your changes
Once an analysis is complete, ensure all your changes are commited and pushed to your project's repository.

Check which files need updating or added to the repository
```sh
git status
```

Stage existing files to be commited, add any addtional files that are needed, then commit those changes including a little description
```sh
git add .
git commit -a -m "Analysis as of commit date"
```

Update the project's repository
```sh
git push -u origin analysis
```

## Finalizing things
Once an analysis is complete, ensure that
- [ ] The delivery folder is uploaded to Box
- [ ] All your changes are committed pushed to your project's repository
