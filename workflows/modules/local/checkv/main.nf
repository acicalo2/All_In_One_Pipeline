#!/usr/bin/ nextflow

nextflow.enable.dsl=2

//  checkv

process CheckV {
    tag {sample_id}
    errorStrategy 'ignore'
    publishDir "${params.outdir}",
                mode: 'copy'

    input:
    tuple val(sample_id),file(contigs_fasta)

    output:
    tuple val(sample_id), file("*"), emit: checkV_ch

    when:
    params.hybrid

    script:
    """
    if [ "${params.dragonflye}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate checkv
        mamba env config vars set CHECKVDB=${params.checkvdb}
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/checkv/dragonflye/
        checkv end_to_end ${params.outdir}/${params.projectID}/${sample_id}/dragonflye/contigs.fa ${params.outdir}/${params.projectID}/${sample_id}/checkv/dragonflye \
               -t ${params.threads}
    elif [ "${params.dragonflye_isolate}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate checkv
        mamba env config vars set CHECKVDB=${params.checkvdb}
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/checkv/dragonflye_isolate/
        checkv end_to_end ${params.outdir}/${params.projectID}/${sample_id}/dragonflye_isolate/contigs.fa ${params.outdir}/${params.projectID}/${sample_id}/checkv/dragonflye_isolate \
               -t ${params.threads}

    elif [ "${params.unicycler}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate checkv
        mamba env config vars set CHECKVDB=${params.checkvdb}
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/checkv/unicycler/
        checkv end_to_end ${params.outdir}/${params.projectID}/${sample_id}/unicycler/assembly.fasta ${params.outdir}/${params.projectID}/${sample_id}/checkv/unicycler \
               -t ${params.threads}
    elif [ "${params.hecatomb}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate checkv
        mamba env config vars set CHECKVDB=${params.checkvdb}
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/checkv/hecatomb/
        checkv end_to_end ${params.outdir}/${params.projectID}/${sample_id}/hecatomb/merged_assembly.fasta ${params.outdir}/${params.projectID}/${sample_id}/checkv/hecatomb \
               -t ${params.threads}
    elif [ "${params.metaspades}" == "true" ];then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate checkv
        mamba env config vars set CHECKVDB=${params.checkvdb}
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/checkv/metaspades/
        checkv end_to_end ${params.outdir}/${params.projectID}/${sample_id}/metaspades/contigs.fasta ${params.outdir}/${params.projectID}/${sample_id}/checkv/metaspades \
               -t ${params.threads}

    elif [ "${params.metaspades}" == "true" ] || [ "${params.phables}" == "true" ] ;then
        eval \"\$(command conda 'shell.bash' 'hook' 2> /dev/null)\"
        conda activate python_mods
    
        mkdir -p ${params.outdir}/${params.projectID}/${sample_id}/post_processing
    
        python3 ${baseDir}/scripts/parse_phables_id_outputs.py \
                -p ${params.outdir}/${params.projectID}/${sample_id}/phables/ \
                -v /${params.outdir}/${params.projectID}/${sample_id}/checkv/ \
                -c ${params.outdir}/${params.projectID}/${sample_id}/metaspades/contigs.fasta \
                -o ${params.outdir}/${params.projectID}/${sample_id}/post_processing/ \
        echo "Process Phage ID Finished" > Process_Phage_ID_Finished 
    
    else 
        "skip this process"
    
    fi
    """
}
