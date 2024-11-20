# Get all analyses
analyses = {v['analysis']: dict(v) for v in config["analyses"]}

rule deseq:
    input:
        samples = lambda wc: analyses[wc.analysis]['samples'],
        counts = countdir + "/{dataset}.csv",
        config = ancient(config["analysis"])
    output:
        xls  = deseqdir + "/{analysis}/{dataset}/analysis.xlsx",
        cnt  = deseqdir + "/{analysis}/{dataset}/counts.csv",
        norm = deseqdir + "/{analysis}/{dataset}/normcounts.csv"
    log:
        deseqdir + "/logs/{analysis}-{dataset}.log"
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb = 8000
    script:
        "../scripts/R/deseq.R"

rule run_deseq:
    input: 
      expand(
        rules.deseq.output,
        analysis=analyses,
        sample=samples,
        dataset=datasets)
