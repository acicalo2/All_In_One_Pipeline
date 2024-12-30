#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*

##################################################################
#			VIRSORTER				 #
##################################################################
*/

process Virsorter {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}",
                mode: 'copy'
    label 'normal'
    input:

    tuple val(sample_id),file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: virsorter_ch

    when:
    params.hybrid
    params.virsorter_run

    script:
    """
    if [ "${params.dragonflye}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate virsorter
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/virsorter/dragonflye
        virsorter  config --set HMMSEARCH_THREADS=${params.virsorter_threads}
        virsorter  config --set CLASSIFY_THREADS=${params.virsorter_threads}
        virsorter run --working-dir ${params.outdir}/${params.projectID}/${sample_id}/virsorter/dragonflye \
                      --seqfile ${params.outdir}/${params.projectID}/${sample_id}/dragonflye/contigs.fa \
                      --verbose --use-conda-off \
                      --db-dir ${params.virsorter_db} 
    elif [ "${params.dragonflye_isolate}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate virsorter
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/virsorter/dragonflye_isolate
        virsorter  config --set HMMSEARCH_THREADS=${params.virsorter_threads}
        virsorter  config --set CLASSIFY_THREADS=${params.virsorter_threads}
        virsorter run --working-dir ${params.outdir}/${params.projectID}/${sample_id}/virsorter/dragonflye_isolate \
                      --seqfile ${params.outdir}/${params.projectID}/${sample_id}/dragonflye_isolate/contigs.fa \
                      --verbose --use-conda-off \
                      --db-dir ${params.virsorter_db} 
    elif [ "${params.unicycler}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate virsorter
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/virsorter/unicycler
        virsorter run --working-dir ${params.outdir}/${params.projectID}/${sample_id}/virsorter/unicycler \
                      --seqfile ${params.outdir}/${params.projectID}/${sample_id}/unicycler/assembly.fasta  \
                      --verbose --use-conda-off \
                      --db-dir ${params.virsorter_db} 
    elif [ "${params.hecatomb}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate virsorter
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/virsorter/hecatomb
        virsorter  config --set HMMSEARCH_THREADS=${params.virsorter_threads}
        virsorter  config --set CLASSIFY_THREADS=${params.virsorter_threads}
        virsorter run --working-dir ${params.outdir}/${params.projectID}/${sample_id}/virsorter/hecatomb \
                      --seqfile ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/*_assembly.fasta  \
                      --verbose --use-conda-off \
                      --db-dir ${params.virsorter_db} \
                      ${params.virsorteropts}
    elif [ "${params.metaspades}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate virsorter
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/virsorter/metaspades
        virsorter run --working-dir virsorter/metaspades \
                      --seqfile ${params.outdir}/${params.projectID}/${sample_id}/metaspades/contigs.fasta  \
                      --verbose --use-conda-off \
                      --db-dir ${params.virsorter_db} \
                      ${params.virsorteropts} \
    else
        "skip this process"
    fi
    """
}