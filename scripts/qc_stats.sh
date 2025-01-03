nextflow run  /export/nextflow/bdrd/all_in_one_pipeline/qc_stats.nf --samplesheet_csv=/mnt/genomics/common_projects/MIDRP/M124007_Thermostable_vaccine_delivery_platform_targeting_members_of_the_Bunyavirales_order/qc_stats/fastq/samplesheet.csv \
                          --outdir=/mnt/genomics/common_projects/MIDRP/M124007_Thermostable_vaccine_delivery_platform_targeting_members_of_the_Bunyavirales_order/qc_stats/fastq/ \
                          --fastqc_threads=100 \
                          --project_id="MI24007" \
                          -profile cluster_normal



#nextflow run  /export/nextflow/bdrd/all_in_one_pipeline/qc_stats.nf --samplesheet_csv=/mnt/genomics/common_projects/NAVSEA/241203_Data_QC/sheet/samplesheet.csv \
#       			  --outdir=/mnt/genomics/common_projects/NAVSEA/241203_Data_QC/pipeline_output/ \
#			  --fastqc_threads=100 \
#			  --project_id="NAVSEA" \
#			  -profile cluster_normal 
