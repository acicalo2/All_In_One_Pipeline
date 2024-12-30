#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  QC 

/* Create fastqc reports */

process Interleave {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/interleave_fastq/", mode: 'copy'
    label 'high'
    //conda './md.yml'

    input:
    tuple val(sample_id), path(fastq_1),path(fastq_2),path(data_dir)

    output:
    tuple val(sample_id), file("*"), emit: interleave_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate bbmap
    mkdir -p ${params.outdir}/${params.project_id}/interleave_fastq/
    reformat.sh in1=${params.fastq_dir}/${fastq_1} in2=${params.fastq_dir}/${fastq_2} out=${sample_id}.fastq.gz ow=t
    """
}

process Pretrim_fastqc_merged {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/pretrim/", mode: 'copy'
    label 'high'
    //conda './md.yml'

    input:
    tuple val(sample_id), file(interleave_files)

    output:

    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md
    mkdir -p ${params.outdir}/${params.project_id}/fastqc/pretrim/
    fastqc --outdir . \
           --threads ${params.fastqc_threads} ${params.outdir}/${params.project_id}/interleave_fastq/*.fastq.gz
    """
}

process Pretrim_fastqc {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/pretrim/", mode: 'copy'
    label 'high'
    //conda './md.yml'

    input: 
    tuple val(sample_id), path(fastq_1),path(fastq_2),path(data_dir)

    output:
    
    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md
    mkdir -p ${params.outdir}/${params.project_id}/fastqc/pretrim/
    fastqc --outdir . \
           --threads ${params.fastqc_threads} ${data_dir}/${sample_id}*.fastq.gz
    """
}

// STAGE 2: Quality control with fastp
process Quality_Control {
    tag {sample_id}
    /* errorStrategy 'ignore' */
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/", mode: 'copy'
    label 'high'
    //conda './vs.yml'
    input: 
    tuple val(sample_id), file(fastq_1), file(fastq_2), file(long_read)

    
    output:
    
    tuple val(sample_id), file("*"), emit: qc_ch

    script:

        """
        if [ "${params.longreads}" == "true" ]; then
            eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
            conda activate vs 
            mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/
            fastp --dedup --overrepresentation_analysis --low_complexity_filter --trim_poly_x \
                  --report_title '${sample_id}_fastp_report' -q 10 -e 10 --cut_front --cut_tail -w 16 \
                  -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \
                  --in1 ${long_read} \
                  -o ${sample_id}.LR.fastp.fastq.gz
        
            echo "Stage 2 Quality Control Finished" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage2_qc.finished 
        
        elif [ "${params.hybrid}" == "true" ]; then 
            eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
            conda activate vs 
            mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/
            fastp --dedup --overrepresentation_analysis --low_complexity_filter --trim_poly_x \
                  --report_title '${sample_id}_fastp_report' -q 10 -e 10 --cut_front --cut_tail -w 16 \
                  -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \
                  --in1 ${long_read} \
                  -o ${sample_id}.LR.fastp.fastq.gz
           
            fastp --verbose -xyp --dedup --report_title "${sample_id}_fastp_report" -q 20 -e 20 --cut_front --cut_tail -w 16 \
                  -j ${params.outdir}/${params.project_id}/${sample_id}/${sample_id}_fastp_ILLUMINA.json \
                  -h ${params.outdir}/${params.project_id}/${sample_id}/${sample_id}_fastp_ILLUMINA.html \
                  -m \
		  --in1 ${fastq_1} --in2 ${fastq_2} \
		  --out1 ${sample_id}_fastp_R1.fastq.gz \
                  --out2 ${sample_id}_fastp_R2.fastq.gz \
                  --merged_out ${sample_id}_fastp_merged.fastq.gz 

            echo "Stage 2 Quality Control Finished" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage2_qc.finished 

        else
            eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
            conda activate vs 
            mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/status_log/
            fastp --verbose -xyp --dedup --report_title "${sample_id}_fastp_report" -q 20 -e 20 --cut_front --cut_tail -w 16 \
		  -j ${params.outdir}/${params.project_id}/${sample_id}/${sample_id}_fastp_ILLUMINA.json \
              	  -h ${params.outdir}/${params.project_id}/${sample_id}/${sample_id}_fastp_ILLUMINA.html \
              	  -m \
		  --in1 ${fastq_1} --in2 ${fastq_2} \
		  --out1 ${sample_id}_fastp_R1.fastq.gz \
              	  --out2 ${sample_id}_fastp_R2.fastq.gz \
              	  --merged_out ${sample_id}_fastp_merged.fastq.gz
            echo "Stage 2 Quality Control Finished" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage2_qc.finished 
        
	    fi
	
	""" 
}

process Post_trim_fastqc {
    tag {sample_id} 
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/fastqc/post_trim/", mode: 'copy'
    label 'high'
    //conda './md.yml'

    input:
    tuple val(sample_id), path(qc_channel_files)

    output:

    tuple val(sample_id), file("*"), emit: posttrim_fastqc_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md
    mkdir -p ${params.outdir}/${params.project_id}/fastqc/post_trim/
    fastqc --outdir . \
           --threads ${params.fastqc_threads} ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/*merged*.fastq.gz

    """
}
// Multiqc - QC Stats

process Multiqc_QC_Stats {
   
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/", mode: 'copy'
    label 'high'
    //conda './md.yml'

    input:
    file fastqc_pre_trim_fles
    file fastqc_post_trim_files 
   
    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate multiqc
    mkdir -p ${params.outdir}/${params.project_id}/multiqc/pretrim/
    multiqc ${params.outdir}/${params.project_id}/fastqc/pretrim/ --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/pretrim/
    mkdir -p ${params.outdir}/${params.project_id}/multiqc/post_trim/
    multiqc ${params.outdir}/${params.project_id}/fastqc/post_trim/ --data-format csv --outdir ${params.outdir}/${params.project_id}/multiqc/post_trim/
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate nextflow
    mkdir -p ${params.outdir}/${params.project_id}/qc_stats/
    Rscript /export/nextflow/bdrd/all_in_one_pipeline/scripts/create_qc_stats.R -i ${params.outdir}/${params.project_id}/multiqc/pretrim/multiqc_data/multiqc_general_stats.csv \
                                                                                -p ${params.outdir}/${params.project_id}/multiqc/post_trim/multiqc_data/multiqc_general_stats.csv \
                                                                                -o ${params.outdir}/${params.project_id}/qc_stats/
    echo "Post Trim Multiqc Complete" > post_trim_multiqc.finished
    """

}

