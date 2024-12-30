#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  hecatomb 

process Hecatomb {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/hecatomb/", mode: 'copy'
    label 'high'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    tuple val(sample_id), path(long_read_fastp)

    output:
    tuple val(sample_id), file("*"), emit: hecatomb_ch

    when:
    params.hecatomb
    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate hecatomb

    mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/hecatomb
    
    hecatomb run all --profile slurm \
    --threads ${threads} \
    --reads ${rpath} \
    --output ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/ \
    ${parmas.hecatomb_config}
    
    echo hecatomb complete > hecatomb_complete.log
    """
}
