rule multiqc:
    input:
        rules.run_fastp.input,
        rules.run_align.input,
    output:
        multiqcdir + "/QC.html"
    log:
        multiqcdir + "/logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    params:
        output_path = multiqcdir,
        input_path = [
          fastqdir + "/json_reports",
          aligndir + "/logs"],
        output_name = "QC.html",
        config = "config/multiqc.yaml"
    shell:
        "(set -x; multiqc -f"
        "  -o {params.output_path} "
        "  -n {params.output_name} "
        "  -c {params.config} "
        "  {params.input_path} "
        ") &> {log}"

rule run_multiqc:
    input: rules.multiqc.output
