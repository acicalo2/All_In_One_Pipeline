#!/usr/bin/ nextflow

nextflow.enable.dsl=2
/*
#########################################################################
                            Move Output
*/

process Move_Output {
	input:
	val input_val

        script:
	"""
	mv -f ${params.outdir}/${params.projectID} ${params.final_data_dir}
	"""
}