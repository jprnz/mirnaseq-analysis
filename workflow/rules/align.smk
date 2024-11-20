# Order of alignments
dataset_order = ["mature", "hairpin"]

# FASTAs used for indexing
datasets = {
  "mature":  aligndir + "/fasta/mature.fasta",
  "hairpin": aligndir + "/fasta/hairpin.fasta",
  # "pirna":   aligndir + "/fasta/pirna.fasta",
  # "ncrna":   ncrna_fasta,
  # "genome":  genome_fasta
}

# Dataset FASTA files (to be filtered)
def get_dataset_fasta(wc):
  if wc.dataset == "mature":
    ret = resource_path / mirna_config['mature']
  elif wc.dataset == "hairpin":
    ret = resource_path / mirna_config['hairpin']
  elif wc.dataset == "pirna":
    ret = resource_path / pirna_config['fasta']
  elif wc.dataset == "ncdna":
    ret = ncdna_fasta
  elif wc.dataset == "genome":
    ret = genome_fasta
  else:
    raise ValueError(f"Dataset not known {wc.dataset}")
  return ret

# Get regex for filtering mi/pirna datasets
def get_dataset_regex(wc):
  if wc.dataset == "mature":
    ret = mirna_dataset['mature']['regex']
  elif wc.dataset == "hairpin":
    ret = mirna_dataset['hairpin']['regex']
  elif wc.dataset == "pirna":
    ret = pirna_dataset['regex']
  else:
    raise ValueError(f"Dataset not known {wc.dataset}")
  return ret

# FASTQ used for alignment
def get_dataset_fastq(wc):
  ind = dataset_order.index(wc.dataset)
  if ind == 0:
    ret = f"{fastqdir}/{wc.sample}.fastq.gz"
  else:
    dat = dataset_order[ind - 1]
    ret = f"{aligndir}/unaligned/{wc.sample}-{dat}.fastq"
  return ret

# Index used for alignment
def get_dataset_index(wc):
  if wc.dataset == "genome":
    ret = genome_index
  elif wc.dataset == "ncrna":
    ret = ncrna_index
  else:
    ret = f"{aligndir}/index/{wc.dataset}/index.1.ebwt"
  return ret

rule make_mirna_fasta:
  input:
    ancient(get_dataset_fasta)
  output:
    aligndir + "/fasta/{dataset}.fasta"
  log:
    aligndir + "/logs/seqkit-{dataset}.log"
  conda:
    "../envs/seqkit.yaml"
  resources:
    mem_mb = 8000
  shell:
    r"(set -x; "
    r"  seqkit grep -i -r -n -p '^{species}(-|_)' {input} "
    r"  | seqkit seq --rna2dna"
    r"  | sed -r 's/>([^\S\|]+).*/>\1/'"
    r"  | seqkit rmdup -n > {output}) &> {log}"

rule bowtie_index:
  input:
    ancient(lambda wc: datasets[wc.dataset])
  output:
    aligndir + "/index/{dataset}/index.1.ebwt"
  log:
    aligndir + "/logs/index-{dataset}.log"
  params:
    path = lambda wc: aligndir + f"/index/{wc.dataset}/index"
  conda:
    "../envs/bowtie.yaml"
  resources:
    mem_mb = 16000
  shell:
    "(set -x; bowtie-build {input} {params.path}) &> {log}"

rule bowtie_align:
  input:
    fq =  ancient(get_dataset_fastq),
    idx = ancient(get_dataset_index)
  output:
    bam = aligndir + "/{sample}-{dataset}.bam",
    unal = aligndir + "/unaligned/{sample}-{dataset}.fastq"
  log:
    aligndir + "/logs/{sample}-{dataset}.log"
  params:
    index = lambda wc, input: input.idx.removesuffix(".1.ebwt"),
    tmpdir = lambda wc: f"{aligndir}/.temp/{wc.dataset}/{wc.sample}"
  conda:
    "../envs/bowtie.yaml"
  resources:
    mem_mb = 8000
  threads: 8
  shell:
    "(set -x;"
    "  [[ -e {params.tmpdir} ]] && rm -r {params.tmpdir};"
    "  mkdir -p {params.tmpdir};"
    "  bowtie"
    "    -v 2"
    "    -k 1"
    "    --best"
    "    --tryhard"
    "    --no-unal"
    "    --chunkmbs 2056"
    "    --mm --time -S"
    "    --threads 32"
    "    --un {output.unal}"
    "    -x {params.index}"
    "    {input.fq}"
    "  | samtools sort"
    "    -u -o '{output.bam}##idx##{output.bam}.bai'"
    "    -T {params.tmpdir}"
    "    --write-index;"
    "  [[ -e {params.tmpdir} ]] && rm -r {params.tmpdir}"
    ") &> {log}"

rule run_align:
  input:
    expand(rules.bowtie_align.output, dataset=datasets, sample=samples)

