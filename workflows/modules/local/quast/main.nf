#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  quast

process Quast {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}", mode: 'copy'
    label 'low'

    input:

    output:
    tuple val(sample_id), file("*"), emit: quast_ch

    script:

    """
    if [ "${params.hecatomb}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate quast

        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/quast

        quast ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
              ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
              ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
             --output-dir ${params.outdir}/${params.projectID}/${sample_id}/quast \
             --threads ${params.threads} \
             --circos ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/merged_assembly.fasta
        
        echo "quast finished" > quast_complete.log
    elif [ "${params.dragonflye}" == "true"]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate quast
        quast --nanopore ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
              --pe1 ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
              --pe2 ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
              --output-dir ${params.outdir}/${params.projectID}/${sample_id}/quast/dragonflye \
              --threads 64
              --circos ${params.outdir}/${params.projectID}/${sample_id}/dragonflye/contigs.fa
        echo "quast finished" > quast_complete.log

    elif [ "${params.unicycler}" == "true"]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate quast

        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/quast/unicycler
        quast --nanopore ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
              --pe1 ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
              --pe2 ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
              --output-dir ${params.outdir}/${params.projectID}/${sample_id}/quast/unicycler \
              --threads 64
              --circos ${params.outdir}/${params.projectID}/${sample_id}/unicycler/contigs.fa
        echo "quast finished" > quast_complete.log

    elif [ "${params.metaspades}" == "true" ] || [ "${params.hybrid}" == "true"];
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate quast

        echo "add the quast script here" > script_needs_added.txt
    else 
      "skip this process"
  
    fi
    """


}
