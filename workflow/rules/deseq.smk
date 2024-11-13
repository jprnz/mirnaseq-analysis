# Get all analyses
analyses = {v['analysis']: dict(v) for v in config["analyses"]}

deseq_targets = list()
for analysis in analyses:
  deseq_targets += expand(
    deseqdir + f"/{analysis}" + "/{dataset}.xlsx",
    dataset = analyses[analysis]['datasets'])

rule deseq:
    input:
        samples = lambda wc: analyses[wc.analysis]['samples'],
        counts = countdir + "/{dataset}.csv",
        config = ancient("config/analysis.yaml"),
    output:
        xls  = deseqdir + "/{analysis}/{dataset}-analysis.xlsx",
        cnt  = deseqdir + "/{analysis}/{dataset}-counts.csv",
        norm = deseqdir + "/{analysis}/{dataset}-normcounts.csv"
    log:
        deseqdir + "/logs/{analysis}-{dataset}.log"
    params:
        analysis = lambda wc: analyses[wc.analysis]
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb = 8000
    script:
        "../scripts/R/deseq.R"

rule run_deseq:
    input: deseq_targets
