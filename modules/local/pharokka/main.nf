#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  pharokka

process Pharokka {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}",
                mode: 'copy'

    input:
    tuple val(sample_id),file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: pharokka_ch

    when:
    params.hybrid

    script:
    """
    if [ "${params.dragonflye}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate pharokka
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/pharokka/dragonflye
        pharokka.py --force --infile ${params.outdir}/${params.projectID}/${sample_id}/dragonflye/contigs.fa \
                    --outdir ${params.outdir}/${params.projectID}/${sample_id}/pharokka/dragonflye \
                    --database ${params.pharokka_db} \
                    --dnaapler \
                    --threads ${threads} \
                    --meta \
                    --meta_hmm
    elif [ "${params.dragonflye_isolate}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate pharokka
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/pharokka/dragonflye_isolate
        pharokka.py --force --infile ${params.outdir}/${params.projectID}/${sample_id}/dragonflye_isolate/contigs.fa \
                    --outdir ${params.outdir}/${params.projectID}/${sample_id}/pharokka/dragonflye_isolate \
                    --database ${params.pharokka_db} \
                    --dnaapler \
                    --threads ${threads} \
                    --meta \
                    --meta_hmm

    elif [ "${params.unicycler}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate pharokka
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/pharokka/unicycler
        pharokka.py --force --infile ${params.outdir}/${params.projectID}/${sample_id}/unicycler/assembly.fasta \
                    --outdir ${params.outdir}/${params.projectID}/${sample_id}/pharokka/unicycler \
                    --database ${params.pharokka_db} \
                    --dnaapler \
                    --threads ${threads} \
                    --meta \
                    --meta_hmm

    elif [ "${params.hecatomb}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate pharokka
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/pharokka/hecatomb
        pharokka.py --force --infile ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/*_assembly.fasta \
                    --outdir ${params.outdir}/${params.projectID}/${sample_id}/pharokka/hecatomb \
                    --database ${params.pharokka_db} \
                    --dnaapler \
                    --threads ${threads} \
                    --meta \
                    --meta_hmm
    elif [ "${params.metaspades}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate pharokka
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/pharokka/metaspades
        pharokka.py --force --infile ${params.outdir}/${params.projectID}/${sample_id}/metaspades/contigs.fasta \
                    --outdir ${params.outdir}/${params.projectID}/${sample_id}/pharokka/metaspades \
                    --database ${params.pharokka_db} \
                    --dnaapler \
                    --threads ${params.threads} \
                    --meta \
                    --meta_hmm
    else
        "skip this process"
    fi
    """
   
}
