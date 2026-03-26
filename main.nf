nextflow.enable.dsl=2

params.input1h = "data/dds_forDeseq_transcrSummarized.RData"
params.input4w = "data/4_weeks/dds_forDeseq.RData"
params.output  = "results/"

process LOAD_FILTER {
    publishDir "${params.output}/filtered", mode: 'copy'

    input:
    path file1h
    path file4w

    output:
    path "dds1h.rds"
    path "dds4w.rds"

    script:
    """
    Rscript bin/load_filter.R ${file1h} ${file4w}
    """
}

process DESEQ_PAIRWISE {
    publishDir "${params.output}/pairwise", mode: 'copy'

    input:
    path dds_obj

    output:
    path "pairwise_results/*.csv"

    script:
    """
    Rscript bin/deseq_pairwise.R ${dds_obj}
    """
}

process DESEQ_GROUP {
    publishDir "${params.output}/groupwise", mode: 'copy'

    input:
    path dds_obj

    output:
    path "group_results/*.csv"

    script:
    """
    Rscript bin/deseq_group.R ${dds_obj}
    """
}

process GO_ENRICHMENT {
    publishDir "${params.output}/go", mode: 'copy'

    input:
    path de_files

    output:
    path "go_results/*.csv"

    script:
    """
    Rscript bin/go_enrichment.R ${de_files}
    """
}

workflow {
    Channel
        .of( params.input1h, params.input4w )
        | LOAD_FILTER

    dds4w = LOAD_FILTER.out[1]

    pair = DESEQ_PAIRWISE(dds4w)
    group = DESEQ_GROUP(dds4w)

    GO_ENRICHMENT(pair.concat(group))
}
``
