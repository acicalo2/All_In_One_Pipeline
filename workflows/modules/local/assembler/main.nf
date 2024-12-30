#!/usr/bin/ nextflow

nextflow.enable.dsl=2
/*
#########################################################################
                            Assesmblies
*/
/* metaSPAdes PE trimmed assembly (dependency STAGE 5) */

process MetaSPAdes_Assembly {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/spades/meta_pe_trim", mode: 'copy'
    label 'normal'
    conda './env/spades.yml'

    input:
    
    tuple val(sample_id), val(stage_3b_files)

    output:

    tuple val(sample_id), file("*"), emit: metaspades_assembly_ch

    when:
    params.metaspades

    script:
    """
    if [ "${params.longreads}" == "true" ]; then
        echo "skip this process" > stage6_skipped.txt
    elif [ "${params.hybrid}" == "true" ]; then 
        echo "skip this process" > stage6_skipped.txt
    else
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate spades   
        ${params.metaSPAdes} -t 128 \
                --12 ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr.fastq.gz \
                -o .

        touch ${params.outdir}/${params.project_id}/${sample_id}/status_log/j6_finished_metaSPAdes_assembly.finished 
    fi
    """
} 

/* Dragonflye meta hybrid assembly */
process Dragonflye_Assembly {
    tag {sample_id}
    errorStrategy 'ignore' 
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dragonflye/hybrid/", mode: 'copy'
    label 'normal'

    input:

    tuple val(sample_id), val(stage_3a_files)
    tuple val(sample_id), val(stage_5_files)
    tuple val(sample_id), val(stage_5b_files)

    output:

    tuple val(sample_id), file("*"), emit: dragonflye_assembly_ch

    when: 

    params.dragonflye 

    script:

    """
    if [ "${params.remove_common_flora}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate dragonflye
        dragonflye --force \
                   --cpus ${params.threads} \
                   --tmpdir /dev/shm \
                   --ram ${params.memory} \
                   --racon 4 \
                   --outdir . \
                   --reads ${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_host_contaminant_rRNA_removed_LR.fastq.gz \
                   --R1 ${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_host_contaminant_rRNA_removed_R1.fastq.gz \
                   --R2 ${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_host_contaminant_rRNA_removed_R2.fastq.gz 

        echo "Stage 6b" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage6b.finished 
    else 
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate bbmap    
        reformat.sh in=${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr.fastq.gz \
                    out1=${sample_id}_host_removed_sr_R1.fastq.gz \
                    out2=${sample_id}_host_removed_sr_R2.fastq.gz
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate dragonflye
        dragonflye --force \
                   --cpus ${params.threads} \
                   --tmpdir /dev/shm \
                   --ram ${params.memory} \
                   --racon 4 \
                   --outdir . \
                   --reads ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_LR.fastq.gz \
                   --R1 ${sample_id}_host_removed_sr_R1.fastq.gz \
                   --R2 ${sample_id}_host_removed_sr_R2.fastq.gz

        echo "Stage 6b" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage6b.finished 
    fi
    """
}

/* Dragonflye isolate long assembly */

process Dragonflye_Isolate_Assembly{
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/dragonflye/LR/", mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), val(stage_3a_files)
    tuple val(sample_id), val(stage_5b_files)

    output:
    tuple val(sample_id), file("*"), emit: dragonflye_isolate_ch

    when: 
    params.dragonflye_isolate 

    script:

    """
    if [ "${params.remove_common_flora}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate dragonflye
        dragonflye --force \
                   --cpus ${params.threads} \
                   --tmpdir /dev/shm \
                   --ram ${params.memory} \
                   --racon 4 \
                   --outdir . \
                   --reads ${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_host_contaminant_rRNA_removed_LR.fastq.gz 
   
        echo "Stage 6c" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage6c.finished 
   

    else
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate dragonflye    
        dragonflye --force \
                   --cpus ${params.threads} \
                   --tmpdir /dev/shm \
                   --ram ${params.memory} \
                   --racon 4 \
                   --outdir . \
                   --reads ${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_host_removed_LR.fastq.gz

        echo "Stage 6c" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage6c.finished 
    
    fi
    """
}

/* Unicycler  */

process Unicycler {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/unicycler/", mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), val(stage_3a_files)
    tuple val(sample_id), val(stage_5b_files)

    output:
    tuple val(sample_id), file("*"), emit: unicycler_ch

    when: 
    params.unicycler 

    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate bbmap
    reformat.sh in=${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr.fastq.gz \ 
                out1=${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr_R1.fastq.gz \ 
                out2=${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr_R2.fastq.gz
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate unicycler    
    unicycler -1 ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr_R1.fastq.gz \
              -2 ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr_R2.fastq.gz \
              --min_fasta_length 1000 -t ${params.cpu} \
              --mode normal \
              --out .
    """
}

