rule multiqc:
    input:
        rules.run_fastp.output,
        rules.run_align.output,
        rules.run_dedup.output,
    output:
        multiqcdir + "/QC.html"
    log:
        multiqcdir + "/logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    params:
        input_path = [fastqdir, aligndir + "/logs", dedupdir],
        output_path = multiqcdir,
        config = "config/multiqc.yaml"
    shell:
        "(set -x; multiqc"
        "  --filename QC"
        "  --force"
        "  --verbose"
        "  --outdir {params.output_path}"
        "  --config {params.config} "
        "  {params.input_path}"
        ") &> {log}"

rule run_multiqc:
    input: rules.multiqc.output
