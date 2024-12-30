#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   MMSEQs
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/

/* mmseq contigs */
process MMSEQ_contigs_against_NT {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'
    input:

    tuple val(sample_id), file(blastx_contigs)

    output:

    tuple val(sample_id),file("${sample_id}_mmseq_contig_against_NT.out") ,emit: mmseqs_spades_contigs_ch
    script:
    """
    if [ "${params.metaspades}" == "true" ]; then
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate vs
      
      mmseqs easy-search ${params.outdir}/${params.project_id}/${sample_id}/spades/meta_pe_trim/contigs.fasta \
      ${params.db_BN} \
      ${sample_id}_mmseq_contig_against_NT.out \
      /tmp \
      --max-accept 25 \
      --threads ${params.threads} \
      -s 7.0 -e 1.0E-8 --search-type 3 \
      --split-memory-limit 120G \
      --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen
    
    elif [ "${params.dragonflye}" == "true" ]; then 
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate vs 
      mmseqs easy-search ${params.outdir}/${params.project_id}/${sample_id}/dragonflye/hybrid/contigs.reoriented.fa \
      ${params.db_BN} \
      ${sample_id}_mmseq_contig_against_NT.out \
      /tmp \
      --max-accept 25 \
      --threads ${params.threads} \
      -s 7.0 -e 1.0E-8 --search-type 3 \
      --split-memory-limit 120G \
      --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen

    elif [ "${params.dragonflye_isolate}" == "true" ]; then 
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate vs 
      mmseqs easy-search ${params.outdir}/${params.project_id}/${sample_id}/dragonflye/LR/contigs.reoriented.fa \
      ${params.db_BN} \
      ${sample_id}_mmseq_contig_against_NT.out \
      /tmp \
      --max-accept 25 \
      --threads ${params.threads} \
      -s 7.0 -e 1.0E-8 --search-type 3 \
      --split-memory-limit 120G \
      --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen
    
    elif [ "${params.unicycler}" == "true" ]; then 
      mmseqs easy-search ${params.outdir}/${params.project_id}/${sample_id}/unicycler/normal/assembly.fasta \
      ${params.db_BN} \
      ${sample_id}_mmseq_contig_against_NT.out \
      /tmp \
      --max-accept 25 \
      --threads ${params.threads} \
      -s 7.0 -e 1.0E-8 --search-type 3 \
      --split-memory-limit 120G \
      --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qlen,tlen
    else 
      "skip this process"    
    fi
    """
}
