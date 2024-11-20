rule dedup:
  input:
    aligndir + "/{sample}-{dataset}.bam"
  output:
    dedupdir + "/{sample}-{dataset}.bam"
  log:
    dedupdir + "/logs/{sample}-{dataset}.log"
  conda:
    "../envs/umi-tools.yaml"
  resources:
    mem_mb = 16000
  shell:
      "(set -x; "
      "  jar=\"$(readlink -f $(which umicollapse)).jar\" &&"
      "  java -Xss1G -Xms8G -Xmx16G -jar $jar bam"
      "    --two-pass --umi-separator _"
      "    -i {input}"
      "    -o {output}) &> {log}"

rule run_dedup:
  input: expand(rules.dedup.output, sample=samples, dataset=datasets)
