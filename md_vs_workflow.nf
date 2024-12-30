#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   Metadetector-VS Nextflow
========================================================================================
   Github   : acicalo2
   Contact  :     
----------------------------------------------------------------------------------------

*/


// Include modules
include { Pretrim_fastqc_merged } from './modules/local/qc/main.nf'
include { Quality_Control } from './modules/local/qc/main.nf'
include { Post_trim_fastqc } from './modules/local/qc/main.nf'
include { Multiqc_QC_Stats } from './modules/local/qc/main.nf'
include { Interleave } from './modules/local/qc/main.nf'
include { Remove_host_seq_from_reads } from './modules/local/mapper/main.nf'
include { MetaSPAdes_Assembly } from './modules/local/assembler/main.nf'
include { BlastX_contigs } from './modules/local/blastx/main.nf'
include { BlastX_reads } from './modules/local/blastx/main.nf'
include { Meganize_BlastX_Contigs } from './modules/local/blastx/main.nf'
include { DAA2INFO_contigs_daa_file } from './modules/local/blastx/main.nf'
include { Meganize_BlastX_Reads } from './modules/local/blastx/main.nf'
include { MMSEQ_contigs_against_NT } from './modules/local/mmseqs/main.nf'


workflow {
    sample_id                    = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> row.sample_id } 
//    sample_w_short_reads_ch      = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.fastq_1),file(row.fastq_2) ) } 
//    sample_w_long_reads_ch       = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.long_read)) } 
    sample_w_short_long_reads_ch = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.fastq_1),file(row.fastq_2),file(row.long_read)) }
//    pretrim_ch                   = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.data_dir)) }        
    pretrim_ch                   = Channel.fromPath(params.samplesheet_csv) splitCsv(header:true) .map { row -> tuple(row.sample_id,file(row.fastq_1),file(row.fastq_2),file(row.data_dir)) }      
    Interleave(pretrim_ch)
    Pretrim_fastqc_merged(Interleave.out.interleave_ch)
    Quality_Control(sample_w_short_long_reads_ch)
    Post_trim_fastqc(Quality_Control.out.qc_ch)
    Multiqc_QC_Stats(Pretrim_fastqc_merged.out.pretrim_fastqc_ch.collect(),Post_trim_fastqc.out.posttrim_fastqc_ch.collect())
    Remove_host_seq_from_reads(Quality_Control.out.qc_ch)
    MetaSPAdes_Assembly(Remove_host_seq_from_reads.out.host_removed_ch)
    BlastX_contigs(MetaSPAdes_Assembly.out.metaspades_assembly_ch)
    BlastX_reads(Remove_host_seq_from_reads.out.host_removed_ch)
    Meganize_BlastX_Contigs(BlastX_contigs.out.blastx_contigs_ch)
    Meganize_BlastX_Reads(BlastX_reads.out.blastx_reads_ch)
    MMSEQ_contigs_against_NT(MetaSPAdes_Assembly.out.metaspades_assembly_ch)
    DAA2INFO_contigs_daa_file(Meganize_BlastX_Contigs.out.meganize_blastx_contigs_ch)
}













