#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*

/*

##################################################################
#			PHABLES					 #
##################################################################

*/

process Phables {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.projectID}/${sample_id}",
                mode: 'copy'

    input:
    tuple val(sample_id),file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: phables_ch

    when:
    params.phables
    script:

    """
    if [ "${params.dragonflye}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate phables
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/phables/dragonflye

        phables run --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp
            --input ${params.outdir}/${params.projectID}/${sample_id}/dragonflye/flye-unpolished.gfa \
            --output phables/dragonflye \
            --threads ${params.threads}
        echo "phables finished" > phables.completed

    elif [ "${params.dragonflye_isolate}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate phables
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/phables/dragonflye_isolate

        phables run --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp
            --input ${params.outdir}/${params.projectID}/${sample_id}/dragonflye_isolate/flye-unpolished.gfa \
            --output phables/dragonflye_isolate \
            --threads ${params.threads}
        echo "phables finished" > phables.completed

    elif [ "${params.unicycler}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate phables
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/phables/unicycler

        phables run --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp
                    --input ${params.outdir}/${params.projectID}/${sample_id}/unicycler/assembly.gfa \
                    --output phables/unicycler \
                    --threads ${params.threads}
        echo "phables finished" > phables.completed
    elif [ "${params.hecatomb}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate phables
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/phables/hecatomb

        phables run --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp
            --input ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/merged_assembly.gfa \
            --output phables/hecatomb \
            --threads ${params.threads}        
        echo "phables finished" > phables.completed

    elif [ "${params.metaspades}" == "true" ] || [ "${params.hybrid}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate phables
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/phables/metaspades

        phables run --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp \
                --input ${params.outdir}/${params.projectID}/${sample_id}/metaspades/assembly_graph_with_scaffolds.gfa \
                --output phables/metaspades \
                --threads ${params.threads}    
        echo "phables finished" > phables.completed
    elif [ "${params.metaspades}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate phables
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/phables/metaspades

        phables run --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp \
            --input ${params.outdir}/${params.projectID}/${sample_id}/metaspades/assembly_graph_after_simplification.gfa \
            --output phables/metaspades \
            --threads ${params.threads}
        echo "phables finished" > phables.completed

    else
        "skip this process"
    """

}