#!/usr/bin/env nextflow

process run_deeptrio {
    publishDir params.outDir, mode: 'copy'

    input:
    path reference
    path reference_index
    val child_id
    path child_bam
    path child_bam_index
    val parent_id_1
    path parent_bam_1
    path parent_bam_1_index
    val parent_id_2
    path parent_bam_2
    path parent_bam_2_index


    output:
    path "*.visual_report.html"
    
    path "${child_id}.deeptrio.vcf.gz", emit: child_vcf
	path "${child_id}.deeptrio.vcf.gz.tbi", emit: child_vcf_index
    
    path "${parent_id_1}.deeptrio.vcf.gz", emit: parent_one_vcf
	path "${parent_id_1}.deeptrio.vcf.gz.tbi", emit: parent_one_vcf_index

    path "${parent_id_2}.deeptrio.vcf.gz", emit: parent_two_vcf
	path "${parent_id_2}.deeptrio.vcf.gz.tbi", emit: parent_two_vcf_index
    
    script:
    """
/opt/deepvariant/bin/deeptrio/run_deeptrio \
  --model_type=WGS \
  --ref=${reference}\
  --reads_child=${child_bam} \
  --reads_parent1=${parent_bam_1} \
  --reads_parent2=${parent_bam_2} \
  --output_vcf_child ${child_id}.deeptrio.vcf.gz \
  --output_vcf_parent1 ${parent_id_1}.deeptrio.vcf.gz \
  --output_vcf_parent2 ${parent_id_2}.deeptrio.vcf.gz \
  --sample_name_child ${child_id} \
  --sample_name_parent1 ${parent_id_1} \
  --sample_name_parent2 ${parent_id_2} \
  --num_shards ${task.cpus}  \
  --output_gvcf_child ${child_id}.deeptrio.g.vcf.gz \
  --output_gvcf_parent1 ${parent_id_1}.deeptrio.g.vcf.gz \
  --output_gvcf_parent2 ${parent_id_2}.deeptrio.g.vcf.gz
"""

}



process hiphase {
    publishDir params.outDir, mode: 'copy'
    
    input:
    val sample 
    path deepvariant_vcf
    path deepvariant_vcf_index
    path aligned_bam
    path aligned_bam_index
    path reference
    path reference_index
    val reference_name

    output:
    path "${sample}.deeptrio.phased.vcf.gz", emit: deeptrio_phased_vcf
    path "${sample}.${reference_name}.hiphase.stats.tsv" 
    path "${sample}.${reference_name}.hiphase.blocks.tsv"
    
    shell:
    """
    set -euo pipefail

	hiphase --version

hiphase --threads ${task.cpus} \
--sample-name ${sample} \
--vcf ${deepvariant_vcf} \
--output-vcf "${sample}.deeptrio.phased.vcf.gz" \
--bam ${aligned_bam} \
--reference ${reference} \
--summary-file ${sample}.${reference_name}.hiphase.stats.tsv \
--blocks-file ${sample}.${reference_name}.hiphase.blocks.tsv \
    """
}

workflow {
    deeptrio_ch = run_deeptrio(params.reference, params.reference_index, 
    params.child_id, params.aligned_child_bam, params.aligned_child_bam_index,
    params.parent_id_1, params.aligned_parent_bam_1, params.aligned_parent_bam_index_1,
    params.parent_id_2, params.aligned_parent_bam_2, params.aligned_parent_bam_index_2)
    sample_ch = Channel.of(params.child_id, params.parent_id_1, params.parent_id_2)
    bam_ch = Channel.of(params.aligned_child_bam, params.aligned_parent_bam_1, params.aligned_parent_bam_2)
    bam_index_ch = Channel.of(params.aligned_child_bam_index, params.aligned_parent_bam_index_1, params.aligned_parent_bam_index_2)
    vcf_ch = deeptrio_ch.child_vcf.concat(deeptrio_ch.parent_one_vcf, deeptrio_ch.parent_two_vcf)
    vcf_ch.view()
    vcf_index_ch = deeptrio_ch.child_vcf_index.concat(deeptrio_ch.parent_one_vcf_index, deeptrio_ch.parent_two_vcf_index)
    vcf_index_ch.view()
    hiphase(sample_ch, vcf_ch, vcf_index_ch, bam_ch, bam_index_ch, params.reference, params.reference_index, params.reference_name)
}

