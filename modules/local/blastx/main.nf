#!/usr/bin/ nextflow

nextflow.enable.dsl=2

/*
========================================================================================
   BlastX
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------

*/


/* Post Assembly Steps */ 
process BlastX_contigs {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'
    input:

    tuple val(sample_id), file(stage_6_files)

    output:

    tuple val(sample_id), file("*"), emit: blastx_contigs_ch

    when:

    script:

    """
    if [ "${params.metaspades}" == "true" ]; then
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md  
    diamond blastx ${params.diamond_args} \
            --threads ${params.threads} \
            --db ${params.diamond_dbdir} \
            --query ${params.outdir}/${params.project_id}/${sample_id}/spades/meta_pe_trim/contigs.fasta \
            --range-culling -F 15 \
            --evalue 1e-5 \
            --outfmt 100 --out ${sample_id}_meta_contigs_blastx.daa
    
    echo "blastx spades contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/blastx_spades_contigs.finished 

    elif [ "${params.dragonflye}" == "true" ]; then
    diamond blastx ${params.diamond_args} \
            --threads ${params.threads} \
            --db ${params.diamond_dbdir} \
            --query ${params.outdir}/${params.project_id}/${sample_id}/dragonflye/hybrid/contigs.reoriented.fa \
            --range-culling -F 15 \
            --evalue 1e-5 \
            --outfmt 100 --out ${sample_id}_dragonflye_contigs_blastx.daa
    
    echo "blastx dragonflye contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/blastx_dragonflye_contigs.finished 
    elif [ "${params.dragonflye_isolate}" == "true" ]; then
      diamond blastx ${params.diamond_args} \
              --threads ${params.threads} \
              --db ${params.diamond_dbdir} \
              --query ${params.outdir}/${params.project_id}/${sample_id}/dragonflye/LR/contigs.reoriented.fa \
              --range-culling -F 15 \
              --evalue 1e-5 \
              --outfmt 100 --out ${sample_id}_dragonflye_isolate_contigs_blastx.daa
    
      echo "blastx dragonflye isolate contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/blastx_dragonflye_isolate_contigs.finished 
    elif [ "${params.unicycler}" == "true" ]; then
      diamond blastx ${params.diamond_args} \
              --threads ${params.threads} \
              --db ${params.diamond_dbdir} \
              --query ${params.outdir}/${params.project_id}/${sample_id}/unicycler/normal/assembly.fasta \
              --range-culling -F 15 \
              --evalue 1e-5 \
              --outfmt 100 --out ${sample_id}_unicycler_contigs_blastx.daa
    
      echo "blastx unicycler contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/unicycler_contigs.finished 
    else 
      "skip this process"
   fi 
    """

}

process BlastX_reads {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'

    input:

    tuple val(sample_id), val(reads)

    output:
    tuple val(sample_id), file("*"), emit: blastx_reads_ch
    when:

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate md   
    diamond blastx ${params.diamond_args} \
                    --threads ${params.threads} \
                    --db ${params.diamond_dbdir} \
                    --query ${params.outdir}/${params.project_id}/${sample_id}/trim/quality_control/reads_mapped/${sample_id}_host_removed_sr.fastq.gz \
                    --evalue 1e-5 \
                    --outfmt 100 \
                    --out ${sample_id}_reads_blastx.daa

    echo "blastx reads" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/blastx_reads.finished 

    """
}

/* Meganize BlastX spades contig output */

process Meganize_BlastX_Contigs {
    tag {sample_id}
    //errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'
    conda './env/md.yml'
    input:

    tuple val(sample_id), file(blastx_contigs)

    output:

    tuple val(sample_id), file("*"), emit: meganize_blastx_contigs_ch
    when:

    script:
    """
    if [ "${params.metaspades}" == "true" ]; then
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate md 
      ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
      ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

      ${params.meganpath}/daa-meganizer \
          --in ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx.daa \
          --mapDB ${params.megandb}/megan-map.db \
          --threads ${params.threads} \
          --topPercent 0.5 \
          --minSupportPercent 0 \
          --lcaAlgorithm longReads \
          --longReads true \
          --verbose   

      paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \
        -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx.daa \
        -c2c Taxonomy |  awk '{print \$1}') \
        <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx.daa \
        -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_contigs_blastx_daa_summary_count.tsv
    
      diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
		   --daa ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx.daa >${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx_diamondview.tsv
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate vs

      python /export/nextflow/bdrd/all_in_one_pipeline/scripts/VS_MD_diamond_parser.py -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx_diamondview.tsv -t blastx
      echo "meganize blastx spades contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/meganize_blastx_spades_contigs.finished 
    elif [ "${params.dragonflye}" == "true" ]; then 
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate md       
      ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
      ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

      ${params.meganpath}/daa-meganizer \
          --in ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_contigs_blastx.daa \
          --mapDB ${params.megandb}/megan-map.db \
          --threads ${params.threads} \
          --topPercent 0.5 \
          --minSupportPercent 0 \
          --lcaAlgorithm longReads \
          --longReads true \
          --verbose   

      paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \
        -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_contigs_blastx.daa \
        -c2c Taxonomy |  awk '{print \$1}') \
        <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_contigs_blastx.daa \
        -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_dragonflye_contigs_blastx_daa_summary_count.tsv
    
      diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
		   --daa ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_contigs_blastx.daa > ${sample_id}_dragonflye_contigs_blastx_diamondview.tsv
      echo "meganize blastx dragonflye spades contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/dragonflye_meganize_blastx_spades_contigs.finished 

    elif [ "${params.dragonflye_isolate}" == "true" ]; then 
      ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
      ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

      ${params.meganpath}/daa-meganizer \
          --in ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa \
          --mapDB ${params.megandb}/megan-map.db \
          --threads ${params.threads} \
          --topPercent 0.5 \
          --minSupportPercent 0 \
          --lcaAlgorithm longReads \
          --longReads true \
          --verbose   
        
      paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \
        -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa \
        -c2c Taxonomy |  awk '{print \$1}') \
        <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa \
        -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_dragonflye_isolate_contigs_blastx_daa_summary_count.tsv
    
      diamond view --daa ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa > ${sample_id}_dragonflye_isolate_contigs_blastx_diamondview.tsv
      echo "meganize blastx dragonflye isolate spades contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/dragonflye_isolate_meganize_blastx_spades_contigs.finished 

    elif [ "${params.unicycler}" == "true" ]; then 
      ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
      ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

      ${params.meganpath}/daa-meganizer \
          --in ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_contigs_blastx.daa \
          --mapDB ${params.megandb}/megan-map.db \
          --threads ${params.threads} \
          --topPercent 0.5 \
          --minSupportPercent 0 \
          --lcaAlgorithm longReads \
          --longReads true \
          --verbose    

      paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_contigs_blastx.daa \
      -c2c Taxonomy |  awk '{print \$1}') \
      <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_contigs_blastx.daa \
      -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}unicycler_contigs_blastx_daa_summary_count.tsv

      diamond view --daa ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_contigs_blastx.daa > ${sample_id}_unicycler_contigs_blastx_diamondview.tsv
      echo "meganize blastx dragonflye spades contigs" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/unicycler_meganize_blastx_spades_contigs.finished 
    
    else 
      "skip this process"
  
  fi
    """
}
// DAA2INFO contigs .daa file

process DAA2INFO_contigs_daa_file {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'

    input:

    tuple val(sample_id), file(blastx_contigs)

    output:

    tuple val(sample_id), file("${sample_id}_c2c.txt"), emit: c2c_txt_file


    script:

    """
    if [ "${params.metaspades}" == "true" ]; then
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_meta_contigs_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u

    elif [ "${params.hybrid}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_contigs_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    elif [ "${params.dragonflye_isolate}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_dragonflye_isolate_contigs_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    elif [ "${params.unicycler}" == "true" ]; then 
      bash ${params.meganpath}/daa2info \
      -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_unicycler_contigs_blastx.daa \
      -o ${sample_id}_c2c.txt \
      -c2c Taxonomy -n -r -u
    else 
      "skip this process"
  
    fi
    """
}

/* Meganize Blastx Reads */

process Meganize_BlastX_Reads {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.project_id}/${sample_id}/blast/", mode: 'copy'
    label 'normal'
    input:
    tuple val(sample_id), file(blastx_reads)

    output:

    tuple val(sample_id), file("${sample_id}*.tsv"), emit: meganize_blastx_reads_ch

    when:

    script:
    """
    if [ "${params.map2reads}" == "true" ]; then
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate md  
      ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
      ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

      ${params.meganpath}/daa-meganizer \
      --in ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa \
      --only Taxonomy \
      --mapDB ${params.megandb}/megan-map.db \
      --threads ${params.threads} \
      --minSupportPercent 0 \
      --topPercent 0.5 \
      --lcaAlgorithm weighted \
      --longReads false \
      --verbose

      paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \
        -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa \
        -c2c Taxonomy | awk '{print \$1}') \
        <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa \
        -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_reads_blastx_daa_summary_count.tsv

      diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
		   --daa ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa > ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx_diamondview.tsv
      
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate vs

      python /export/nextflow/bdrd/all_in_one_pipeline/scripts/VS_MD_diamond_parser.py -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx_diamondview.tsv -t blastx

    echo "meganize blastx reads" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/meganize_blastx_reads.finished 
    elif [ "${params.remove_host_seq_from_reads}" == "true" ]; then 
      eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
      conda activate md        
      ln -sf ${params.megandb}/ncbi.map ./ncbi.map 
      ln -sf ${params.megandb}/ncbi.tre ./ncbi.tre 

      ${params.meganpath}/daa-meganizer \
      --in ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa \
      --only Taxonomy \
      --mapDB ${params.megandb}/megan-map.db \
      --threads ${params.threads} \
      --minSupportPercent 0 \
      --topPercent 0.5 \
      --lcaAlgorithm weighted \
      --longReads false \
      --verbose

      paste <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def \
        -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa \
        -c2c Taxonomy | awk '{print \$1}') \
        <(${params.meganpath}/daa2info -P ${params.meganpath}/.MEGAN.def -i ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa \
        -p -c2c Taxonomy | awk '{print \$1,\$2}' FS='\\t' OFS='\\t') > ${sample_id}_reads_blastx_daa_summary_count.tsv

      diamond view --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qlen slen \
		   --daa ${params.outdir}/${params.project_id}/${sample_id}/blast/${sample_id}_reads_blastx.daa > ${sample_id}_reads_blastx_diamondview.tsv

      echo "meganize blastx reads" > ${params.outdir}/${params.project_id}/${sample_id}/status_log/meganize_blastx_reads.finished 
    else 
      "skip this process"
      
    fi
    """  
}
