rule dedup:
  input:
    bam = aligndir + "/{sample}-{dataset}.bam"
  output:
    bam = dedupdir + "/{sample}-{dataset}.bam",
  log:
    dedupdir + "/logs/{sample}-{dataset}.log"
  conda:
    "../envs/umi-tools.yaml"
  resources:
    mem_mb = 16000
  shell:
      "(set -x; umicollapse bam"
      "  --umi-separator _"
      "  -i {input}"
      "  -o {output}) &> {log}"

rule run_dedup:
  input:
    expand(dedupdir + "/{sample}.bam", sample=samples)
