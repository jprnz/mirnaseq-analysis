rule count_genome_saf:
  input:
    dedupdir + "/{sample}-genome.bam"
  output:
    countdir + "/genome/{sample}.saf"
  log:
    countdir + "/logs/genome-saf-{sample}.log"
  conda:
    "../envs/count.yaml"
  shell:
    "(set -x; "
    "  echo -e \"GeneID\tChr\tStart\tEnd\tStrand\" > {output};"
    "  samtools depth {input}"
    "    | awk -v OFS='\t' '{{print $1, $2-1, $2}}'"
    "    | bedtools merge"
    "    | awk -F'\t' -v OFS='\t' '{{"
    "        name = $1 \":\" $2 + 1 \"-\" $3;"
    "        print name, $1, $2 + 1, $3, \"+\""
    "      }}' >> {output}"
    ") &> {log}"

rule count_genome:
  input:
    saf = countdir + "/genome/{sample}.saf",
    bam = dedupdir + "/{sample}-genome.bam"
  output:
    countdir + "/genome/{sample}.tsv"
  log:
    countdir + "/logs/genome-count-{sample}.log"
  params:
    qual = 30
  conda:
    "../envs/count.yaml"
  threads: 8
  shell:
    "(set -x; featureCounts"
    "  -F SAF -d 0 -T 16 -s 0"
    "  -a {input.saf}"
    "  -o {output}"
    "  {input.bam}) &> {log}"

rule count_mirna:
  input:
    dedupdir + "/{sample}-{dataset}.bam"
  output:
    countdir + "/{dataset}/{sample}.tsv"
  log:
    countdir + "/logs/{sample}-{dataset}.log"
  conda:
    "../envs/count.yaml"
  wildcard_constraints:
    dataset = "|".join([v for v in datasets if v != "genome"])
  shell:
    "(set -x; samtools idxstat {input}"
    "  | awk -v OFS='\t' '{{print $1, $3}}' > {output}"
    ") &> {log}"

rule count_gather:
  input:
    expand(countdir + "/{{dataset}}/{sample}.tsv", sample=samples)
  output:
    countdir + "/{dataset}.csv"
  log:
    countdir + "/logs/{dataset}-gather.log"
  script:
    "../scripts/gather_counts.py"

rule run_count:
  input:
    expand(countdir + "/{dataset}.csv", dataset=datasets)

