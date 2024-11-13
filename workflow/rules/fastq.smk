adaptor_5 = config['library']['adaptor_5']
adaptor_3 = config['library']['adaptor_3']
min_read = str(config['library']['min_read'])
umi_length = str(config['library']['umi_length'])

# Default
adaptor_5 = ".*" if not adaptor_5 else adaptor_5
adaptor_3 = ".*" if not adaptor_3 else adaptor_3

bc_pattern = (
    ".+(?P<discard_1>" + adaptor_5 + "){s<=2}"
    "(?P<umi_1>.{" + umi_length + "})"
    "(?P<discard_2>" + adaptor_3 + "){s<=2}.*")

def get_fastqs(wc):
    paths = files.loc[wc.sample].path
    try:
        ret = sorted(paths.tolist())
    except:
        ret = [paths]
    return ret

rule fastq_combine:
    input:
        ancient(get_fastqs)
    output:
        temp(fastqdir + "/fastqs/{sample}.fastq.gz")
    log:
        fastqdir + "/logs/combine-{sample}.log"
    resources:
        mem_mb = 1000
    group: "fastq"
    run:
        if len(input) > 1:
            cmd = "cat {input} > {output}"
        else:
            cmd = "ln -v {input} {output}"
        shell("(set -x; " + cmd + ") &>> {log}")

rule umi_extract:
    input:
        ancient(fastqdir + "/fastqs/{sample}.fastq.gz")
    output:
        temp(fastqdir + "/umi/{sample}.fastq.gz")
    log:
        fastqdir + "/logs/umi-extract-{sample}.log"
    conda:
        "../envs/umi-tools.yaml"
    resources:
        mem_mb = 32000
    group: "fastq"
    shell:
        "(set -x; umi_tools extract"
        "  --stdin {input} "
        "  --stdout {output} "
        "  --extract-method regex "
        "  --bc-pattern '{bc_pattern}' "
        ") &> {log}"

rule fastp:
    input:
        ancient(fastqdir + "/umi/{sample}.fastq.gz"),
    output:
        read = temp(fastqdir + "/{sample}.fastq.gz"),
        json_report = fastqdir + "/json_reports/{sample}.json",
        html_report = fastqdir + "/html_reports/{sample}.html"
    log:
        fastqdir + "/logs/{sample}.log"
    conda:
        "../envs/fastp.yaml"
    resources:
        mem_mb = 16000
    threads: 16
    shell:
        "(set -x; fastp "
        "  -i {input} "
        "  -o {output.read} "
        "  --length_required {min_read} "
        "  --json {output.json_report} "
        "  --html {output.html_report} "
        "  --thread {threads} "
        ") &> {log}"

rule run_fastq_combine:
  input:
      expand(fastqdir + "/fastqs/{sample}.fastq.gz", sample=samples)

rule run_umi_extract:
  input:
      expand(fastqdir + "/umi/{sample}.fastq.gz", sample=samples)

rule run_fastp:
    input:
        expand(fastqdir + "/json_reports/{sample}.json", sample=samples)


