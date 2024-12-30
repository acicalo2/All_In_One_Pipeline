#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Map 2 Reference
========================================================================================

*/

/* map reads to reference sequence and use for downstream analysis (Optional)*/

process Map_Reads_2_RefSeq {
    tag {sample_id}
    errorStrategy 'ignore' 
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/", mode: 'copy'
    label 'normal'
    input: 
    tuple val(sample_id), file(qc_files)

    output:

    tuple val(sample_id), file("*"), emit: map_reads_2_refseq_ch

    when: 
    params.map2reads
    script:
    """
    if [ "${params.longreads}" == "true" ] || [ "${params.targeted}"]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}.LR.fastp.fastq.gz > ${sample_id}_host_removed_LR.sam 
        samtools fastq -f 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_LR.fastq
	
	    pigz *.fastq

        echo "Stage 3a " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage3a.finished 

    elif [ "${params.hybrid}" == "true" ] || [ "${params.targeted}"]; then 
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}.LR.fastp.fastq.gz > ${sample_id}_host_removed_LR.sam 
        samtools fastq -f 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_LR.fastq

        minimap2 -ax sr -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}_fastp_merged.fastq.gz \
        > ${sample_id}_host_removed_sr.sam 

        samtools fastq -f 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_removed_sr.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_sr.fastq
        pigz *.fastq

        echo "Stage 3a " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage3a.finished 
    
    elif [ "${params.targeted}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md 
        minimap2 -ax sr -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}_fastp_merged.fastq.gz \
        > ${sample_id}_host_removed_sr.sam 

        samtools fastq -f 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_removed_sr.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_sr.fastq
        pigz *.fastq

        echo "Stage 3a " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage3a.finished 
    else 
        echo "Stage 3a skipped" > stage3a.skipped
    
    fi
    """

}

/* remove provided host sequence from reads (dependency STAGE 2) */

process Remove_host_seq_from_reads {
    tag {sample_id}
    errorStrategy 'ignore' 
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/", mode: 'copy'
    label 'normal'
    input: 

    tuple val(sample_id), file(mapped_reads_files)

    output:
    
    tuple val(sample_id), file("*"), emit: host_removed_ch

    when: 
    params.remove_host_seq_from_reads
    script:
    if (params.longreads == true)
    """
    if [ "${params.longreads}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}.LR.fastp.fastq.gz > ${sample_id}_host_removed_LR.sam 
        samtools fastq -f 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_LR.fastq
	    pigz *.fastq
        echo "Stage 3b " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage3b.finished 

    elif [ "${params.hybrid}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md 
        minimap2 -ax map-ont -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}.LR.fastp.fastq.gz > ${sample_id}_host_removed_LR.sam 
        samtools fastq -f 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_LR.sam > ${sample_id}_host_LR.fastq

        minimap2 -ax sr -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}_fastp_merged.fastq.gz \
        > ${sample_id}_host_removed_sr.sam 

        samtools fastq -f 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_removed_sr.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_sr.fastq
        pigz *.fastq

        echo "Stage 3b " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage3b.finished 


    else
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md  
        minimap2 -ax sr -t ${params.threads} ${params.host_db_human} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/${sample_id}_fastp_merged.fastq.gz \
        > ${sample_id}_host_removed_sr.sam 

        samtools fastq -f 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_removed_sr.fastq
        samtools fastq -F 4 ${sample_id}_host_removed_sr.sam > ${sample_id}_host_sr.fastq
        pigz *.fastq

        echo "Stage 3b " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage3b.finished 
    fi
    """
}

/* STAGE 4: remove common lab contaminant reads  */


process Remove_Common_Lab_Contaminants {
    tag {sample_id}
    errorStrategy 'ignore' 
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/", mode: 'copy'
    label 'normal'

    input: 
    tuple val(sample_id), file(mapped_reads_files)
    tuple val(sample_id), file(contaminant_removed_files)


    output:
    tuple val(sample_id), file("*"), emit: contaminant_removed_ch

    when: 
    
    params.remove_common_flora 

    script:
    """
    if [ "${params.longreads}" == "true" ] || [ "${params.targeted}"]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -N 1 -k 14 -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_LR.fastq.gz \
        > ${sample_id}_contaminant_removed_LR.sam         
        
        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_LR.sam  > ${sample_id}_contaminant_removed_LR.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_LR.sam > ${sample_id}_contaminant_LR.fastq

        pigz *.fastq

        echo "Stage 4 " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage4.finished 

    elif [ "${params.hybrid}" == "true" ] || [ "${params.targeted}"]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -N 1 -k 14 -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_LR.fastq.gz \
        > ${sample_id}_contaminant_removed_LR.sam   

        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_LR.sam  > ${sample_id}_contaminant_removed_LR.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_LR.sam > ${sample_id}_contaminant_LR.fastq

        minimap2 -ax sr -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_sr.fastq.gz \
        > ${sample_id}_contaminant_removed_sr.sam 

        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_sr.sam  > ${sample_id}_contaminant_removed_sr.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_sr.sam > ${sample_id}_contaminant_sr.fastq
        
        pigz *.fastq

        echo "Stage 4 " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage4.finished 

    elif [ "${params.longreads}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_removed_LR.fastq.gz \
        > ${sample_id}_host_contaminant_removed_LR.sam

        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_LR.sam  > ${sample_id}_contaminant_removed_LR.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_LR.sam > ${sample_id}_contaminant_LR.fastq

        pigz *.fastq

        echo "Stage 4 " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage4.finished 

    elif [ "${params.hybrid}" == "true" ]; then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_removed_LR.fastq.gz \
        > ${sample_id}_host_contaminant_removed_LR.sam   

        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_LR.sam  > ${sample_id}_contaminant_removed_LR.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_LR.sam > ${sample_id}_contaminant_LR.fastq

        minimap2 -ax sr -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_removed_sr.fastq.gz \
        > ${sample_id}_host_contaminant_removed_sr.sam 
    
        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_sr.sam  > ${sample_id}_contaminant_removed_sr.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_sr.sam > ${sample_id}_contaminant_sr.fastq

        echo "Stage 4 " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage4.finished 

    else
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md 
        minimap2 -ax sr -t ${params.threads} ${params.bbmap_ref_contam} \
        ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/targeted_read_mapping/${sample_id}_host_removed_sr.fastq.gz \
        > ${sample_id}_host_contaminant_removed_sr.sam 

        samtools fastq -f 4 \
        ${sample_id}_contaminant_removed_sr.sam  > ${sample_id}_contaminant_removed_sr.fastq
        
        samtools fastq -F 4 \
        ${sample_id}_contaminant_removed_sr.sam > ${sample_id}_contaminant_sr.fastq

        pigz *.fastq

        echo "Stage 4 " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage4.finished 

    fi
    """
    
}

/* remove common flora rRNA reads (dependency STAGE 4) */

process Remove_Common_Flora_rRNA_reads {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/", mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), file(contaminant_removed_files)

    output:
    tuple val(sample_id), file("*"), emit: common_flora_removed_ch

    when:

    params.remove_common_flora 

    script:
    """
    if [ "${params.longreads}" == "true" ]; then
        echo "skip this process" > stage5_processed_skipped.txt
    else
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        bbmap.sh ${params.bbmap_args} \
                 in=${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_contaminant_removed_sr.fastq.gz \
                 path=${params.bbmap_ref_silva} tossbrokenreads printunmappedcount=t \
                 covstats=${sample_id}.rRNA_covstats.txt \
                 outm=${sample_id}_rRNA.fastq.gz usejni=t \
                 -Xmx${params.memory}g \
                 outu=${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz overwrite=true

                 seqtk seq -1 ${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz | /bin/gzip -1 > ${sample_id}_host_contaminant_rRNA_removed_R1.fastq.gz
                 seqtk seq -2 ${sample_id}_host_contaminant_rRNA_removed_pe.fastq.gz | /bin/gzip -1 > ${sample_id}_host_contaminant_rRNA_removed_R2.fastq.gz
    
        echo "Stage 5 " > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage5.finished 
    
    fi
    """
}

/* remove common flora rRNA reads (dependency STAGE 4) */

process Remove_Common_Flora_rRNA_reads_b {
    tag {sample_id}
    errorStrategy 'ignore' 
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/", mode: 'copy'
    label 'normal'
    input:

    tuple val(sample_id), file(contaminant_removed_files)

    output:
    tuple val(sample_id),file("*"), emit: remove_common_flora_rRNA_reads_b_ch

    when: 
    params.hybrid | params.longreads
    script:
    
    """
    if [ "${params.remove_common_flora} == "true"); then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate md
        minimap2 -ax map-ont -t ${params.threads} ${params.host_db_silva} ${params.outdir}/${params.project_id}/${sample_id}/mapped_reads/contaminant_removed/${sample_id}_host_contaminant_removed_LR.fastq > ${sample_id}_host_contaminant_rRNA_removed_LR.sam
        samtools fastq -f 4 ${sample_id}_host_contaminant_rRNA_removed_LR.sam > ${sample_id}_host_contaminant_rRNA_removed_LR.fastq
        samtools fastq -F 4 ${sample_id}_host_contaminant_rRNA_removed_LR.sam > ${sample_id}_rRNA_LR.fastq

        echo "Stage 5b" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage5b.finished 
    else 
        echo "skip this process" > stage5b_processed_skipped.txt
    fi
    """
}
