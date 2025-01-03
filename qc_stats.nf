#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   qc stats Nextflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

// Include modules
include { Pretrim_fastqc_merged } from './modules/local/qc/main.nf'
include { Quality_Control } from './modules/local/qc/main.nf'
include { Post_trim_fastqc } from './modules/local/qc/main.nf'
include { Multiqc_QC_Stats } from './modules/local/qc/main.nf'
include { Interleave } from './modules/local/qc/main.nf'
workflow {
    sample_id                    = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> row.sample_id } 
//    sample_w_short_reads_ch      = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.fastq_1),file(row.fastq_2) ) } 
//    sample_w_long_reads_ch       = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.long_read)) } 
    
    sample_w_short_long_reads_ch = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.fastq_1),file(row.fastq_2),file(row.long_read),
														      file(row.data_dir)) }    
    sample_ch                    = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,row.fastq_1,row.fastq_2,
															      	    row.long_read) }     
    //data_dir                   = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> file(row.data_dir) }
    Interleave(sample_ch)
    Pretrim_fastqc_merged(sample_ch,Interleave.out.interleave_ch)
    Quality_Control(sample_ch)
    Post_trim_fastqc(Quality_Control.out.qc_ch)
    Multiqc_QC_Stats(Pretrim_fastqc_merged.out.pretrim_fastqc_ch.collect(),Post_trim_fastqc.out.posttrim_fastqc_ch.collect())
}
