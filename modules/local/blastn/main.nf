#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  Blastn

process Blastn{
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.projectID}/${sample_id}/blast/",
                mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), path(parse_blastx_output_file)

    output:
    tuple val(sample_id),file("*"), emit: blastn_ch

    script: 
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/blast/

    blastn -query ${params.outdir}/${params.projectID}/${sample_id}/blast/${sample_id}_phables_noIntegrase.fasta \
           -db ${params.blast_phagedb} \
           -task megablast \
           -evalue 1e-8 \
           -max_target_seqs 10 \
           -outfmt 0 \
	   -sorthits 2 \
           -num_threads=${params.threads} \
           -out ${sample_id}_phables_noIntegrase_phages_megablast.out
    """
}

process Parse_Blast_Output{
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}/${params.projectID}/${sample_id}/blast/",
                mode: 'copy'
    label 'normal'

    input:
    tuple val(sample_id), path(process_blastn_file)

    output:
    tuple val(sample_id),file("*"), emit: parse_blastn_output_ch

    script:
    """
    eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
    conda activate vs

    python3 ${baseDir}/scripts/parse_megablast_txt.py \
        -i ${params.outdir}/${params.projectID}/${sample_id}/blast/${sample_id}_phables_noIntegrase_phages_megablast.out \
        -q ${params.outdir}/${params.projectID}/${sample_id}/blast/${sample_id}_phables_noIntegrase.fasta \
        -o ${sample_id}_phables_noIntegrase_final_phages.fasta 
    """
}