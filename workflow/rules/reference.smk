rule download_sequence:
  output:
    f"{ref_path}/{ref_file}" + ".{reference}.fa"
  params:
    **ref_args,
    datatype = "{reference}" 
  log:
    ref_log + "/download-{reference}.log"
  wrapper:
    "v4.7.3/bio/reference/ensembl-sequence"

rule download_annotation:
  output:
    f"{ref_path}/{ref_file}" + ".gtf"
  params:
    **ref_args
  log:
    ref_log + "/download-gtf.log"
  wrapper:
    "v4.7.3/bio/reference/ensembl-annotation"

rule bowtie_index_reference: 
  input:
    f"{ref_path}/{ref_file}" + ".{reference}.fa"
  output:
    ref_path + "/bowtie2-index/{reference}.1.bt2"
  log:
    ref_log + "/bowtie2-index-{reference}.log"
  params:
    path = lambda wc, output: output[0].removesuffix(".1.bt2")
  resources:
    mem_mb = 16000
  threads: 16
  conda:
    "../envs/bowtie.yaml"
  shell:
    "(set -x; "
    "  outpath=$(dirname {output});"
    "  [[ -e $outpath ]] || mkdir -p $outpath;"
    "  bowtie2-build --threads {threads} {input} {params.path}"
    ") &> {log}"

rule make_genes_file:
  input:
    genome_gtf
  output:
    f"{genome_gtf}.genes"
  log:
    ref_log + "/make_genes_file.log"
  script:
    "../scripts/gene-annotations.py"

rule make_faidx:
  input:
    ref_path + "/{reference}.fa"
  output:
    ref_path + "/{reference}.fa.fai"
  log:
    ref_log + "/make_faidx-{reference}.fa"
  conda:
    "../envs/samtools.yaml"
  shell:
    "(set -x; samtools faidx {input}) &> {log}"

