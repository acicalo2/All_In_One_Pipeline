#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  cenotetaker3

process Cenotetaker3 {
    tag {sample_id}
    errorStrategy 'ignore'   
    publishDir "${params.outdir}",
                mode: 'copy'

    input:
    
    tuple val(sample_id),file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: ct3_ch

    when:
    params.hybrid

    script 
    """
    if [ "${params.dragonflye}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate cenote-taker3
        mamba env config vars set CENOTE_DBS=${params.cenote_taker3_db}
        echo "ct3_`date +%Y-%m-%d_%H%M%S`" > date.txt
        timestamp="\$(cat date.txt)"
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/dragonflye
        cenotetaker3 --run_title ct3_\${timestamp} \
                     --contigs ${params.outdir}/${params.projectID}/${sample_id}/dragonflye/contigs.fa \
                     --working_directory ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/dragonflye \
                     --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
                     --prune_prophage F \
                     --hhsuite_tool hhblits \
                     --cpu ${params.threads} \
                    --caller adaptive
        echo "centotaker3 finished" > centotaker3_completed.log

    elif ["${params.dragonflye_isolate}" == "true"];then 
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate cenote-taker3
        mamba env config vars set CENOTE_DBS=${params.cenote_taker3_db}
        echo "ct3_`date +%Y-%m-%d_%H%M%S`" > date.txt
        timestamp="\$(cat date.txt)"
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/dragonflye_isolate
        cenotetaker3 --run_title ct3_\${timestamp} \
                     --contigs ${params.outdir}/${params.projectID}/${sample_id}/dragonflye_isolate/contigs.fa \
                     --working_directory ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/dragonflye_isolate \
                     --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
                     --prune_prophage F \
                     --hhsuite_tool hhblits \
                     --cpu ${params.threads} \
                     --caller adaptive
        echo "centotaker3 finished" > centotaker3_completed.log

    elif ["${params.unicycler}" == "true"];then 
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate cenote-taker3
        mamba env config vars set CENOTE_DBS=${params.cenote_taker3_db}
        echo "ct3_`date +%Y-%m-%d_%H%M%S`" > date.txt
        timestamp="\$(cat date.txt)"
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/unicycler
        cenotetaker3 --run_title ct3_\${timestamp} \
                     --contigs ${params.outdir}/${params.projectID}/${sample_id}/unicycler/assembly.fasta  \
                     --working_directory ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/unicycler \
                     --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
                     --prune_prophage F \
                     --hhsuite_tool hhblits \
                     --cpu ${params.threads} \
                     --caller adaptive
        echo "centotaker3 finished" > centotaker3_completed.log
    
    elif ["${params.hecatomb}" == "true"];then 
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate cenote-taker3
        mamba env config vars set CENOTE_DBS=${params.cenote_taker3_db}
        echo "ct3_`date +%Y-%m-%d_%H%M%S`" > date.txt
        timestamp="\$(cat date.txt)"
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/hecatomb
        cenotetaker3 --run_title ct3_\${timestamp} \
                     --contigs ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/merged_assembly.fasta \
                     --working_directory ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/hecatomb \
                     --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
                     --prune_prophage F \
                     --hhsuite_tool hhblits \
                     --cpu ${params.threads} \
                     --caller adaptive
        echo "centotaker3 finished" > centotaker3_completed.log

    elif ["${params.metaspades}" == "true"];then 
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate cenote-taker3
        mamba env config vars set CENOTE_DBS=${params.cenote_taker3_db}
        echo "ct3_`date +%Y-%m-%d_%H%M%S`" > date.txt
        timestamp="\$(cat date.txt)"
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/metaspades
        cenotetaker3 --run_title ct3_\${timestamp} \
                     --contigs ${params.outdir}/${params.projectID}/${sample_id}/metaspades/contigs.fasta \
                     --working_directory ${params.outdir}/${params.projectID}/${sample_id}/cenotetaker3/metaspades \
                     --reads ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_ONT.fastq.gz ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R1.fastq.gz \
                             ${params.outdir}/${params.projectID}/${sample_id}/fastp/${sample_id}_fastp_R2.fastq.gz \
                     --prune_prophage F \
                     --hhsuite_tool hhblits \
                     --cpu ${params.threads} \
                     --caller adaptive        
        echo "centotaker3 finished" > centotaker3_completed.log
    else 
      "skip this process"
      
    fi
    """
}