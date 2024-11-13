# General settings
To configure this workflow, modify `config/config.yaml` according to your needs paying special attention to the reference genome, output folders, and strandedness options.

# Mapping and QC
FASTQ files should be listed in `fastqs.tsv` file. For each sample, the columns `sample`, `pair`, and `path` must be defined. `sample` indicates a unique unit that will be processed, `pair` indicates the sequence pair of the FASTQ file, and `path` should describe the location of the FASTQ file (as an absolute path or relative to the project directory). Multiple values for the same `sample` and `pair` will merged together before processing,useful when a samples need to be merged together.

Below is an example of how this file might be compiled from the project directory, where `data` contains FASTQ files:
```sh
find data/ -name "*.fastq.gz" \
  | sed -r 's|.*/(.*)_S[0-9]+_.*_(R[12])_.*.fastq.gz|\1\t\2\t\0|g' \
  | sort -V >> config/fastqs.tsv
```

### Strandedness
The `strandedness` option in `config.yaml` defines the strandedness for all samples in the project. However, per-sample strandedness can be defined in a separate file. To do so, include a tab-separated table with the columns `sample` and `strandedness`. Values for this table can be `none`, `yes`, or `reverse` (see `config.yaml` for a description of each of these).


# Differential expression
## Samples and conditions
Samples to be included in an analysis need to be defined in a tab-separated file that contain `sample` (as defined in `fastqs.tsv`) and at least one condition to test. Any combination of samples can be included. An optional column `alias` can be included to change how these names appear downstream of DESeq. Additional columns will be converted to factors and matching the level ordering to what was givein in the input. When doing an all-vs-all comparison using `[<condtion>, all]` as a value for `contrasts`, the sample order will determine which comparisons will be made. Generally it is best to include WT / baseline samples first, followed by the most interesting conditions.

## Analysis configuration
Options to DESeq to regarding experiment design and which conditions to test are setup via the analysis.yaml file:
```yaml
analysese:
  - analysis: name of the analysis
    design: ~ analysis formula
    samples: file containing sample / condition information (ex. samples.tsv)
    group_counts: variable used to calculate mean counts in results table
    pca:
        pcs: list of two PCs to plot
        color_by: variable used to determine color
        shape_by: variable used to determine shape
        label_by: variable used to determine label
        color_legend: name to use for color legend
        shape_legend: name to use for shape legend
        label_legend: name to use for label legend
        alpha: transparency of points, as a percent
        size: size of points
    heatmap:
        only_contrasts: Only use samples defined in contrasts (true or false, default false)
    deseq:
        parameter: value
    results:
      contrasts:
        - [condition, all]
        - <display_name>:
          # Single strings are allowed, but will be coerced to a list
          - [list of coefficients]
          - [list of coefficients]
        - <display_name>:
          - [list of coefficents]

  # Additional analyses can be included... 
  - analysis: Name of second analysis
    samplesheet: subset.tsv
    design: ~ condition
    etc...
```

### Notes
* `results`, and `deseq` can be used to pass any arbitrary parameter to either `results()` or `DESeq()`
* Each element of `contrasts` will be passed to the `results()` function, see the [DESeq2 manual](https://bioconductor.org/packages/devel/bioc/manuals/DESeq2/man/DESeq2.pdf#Rfn.results.1).
* Output of resultsNames() is printed in the log file for the analysis.
* For LRT, `design` is the full model, use the `deseq` to set the reduced model and make sure `test` is `LRT`, and set the value for `contrast` explicitly:
```yaml
design: ~ batch + condition
deseq:
  reduced: ~ batch
  test: LRT
contrasts:
  - LRT:
    - condition
```

