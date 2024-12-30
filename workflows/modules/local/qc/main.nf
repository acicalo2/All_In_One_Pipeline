#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  QC 

/* Create fastqc reports */

process Pretrim_fastqc {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/fastqc/pretrim/", mode: 'copy'
    label 'low'
    //conda './md.yml'

    input: 
    tuple val(sample_id), path(data_dir)

    output:
    
    tuple val(sample_id), file("*"), emit: pretrim_fastqc_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md
    mkdir -p ${params.outdir}/${params.project_id}/${sample_id}/fastqc/pretrim/
    fastqc --outdir . \
           --threads 32 ${data_dir}/[!Undetermined_]*.fastq.gz
    """
}

// STAGE 2: Quality control with fastp
process Quality_Control {
    tag {sample_id}
    /* errorStrategy 'ignore' */
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/", mode: 'copy'
    label 'normal'
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
                  -j ${params.outdir}/${params.projectID}/${sample_id}/${sample_id}_fastp_ILLUMINA.json \
                  -h ${params.outdir}/${params.projectID}/${sample_id}/${sample_id}_fastp_ILLUMINA.html \
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
		      -j ${params.outdir}/${params.projectID}/${sample_id}/${sample_id}_fastp_ILLUMINA.json \
              -h ${params.outdir}/${params.projectID}/${sample_id}/${sample_id}_fastp_ILLUMINA.html \
              -m \
		      --in1 ${fastq_1} --in2 ${fastq_2} \
		      --out1 ${sample_id}_fastp_R1.fastq.gz \
              --out2 ${sample_id}_fastp_R2.fastq.gz \
              --merged_out ${sample_id}_fastp_merged.fastq.gz

            echo "Stage 2 Quality Control Finished" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/stage2_qc.finished 
        
	    fi
	
	""" 
}
