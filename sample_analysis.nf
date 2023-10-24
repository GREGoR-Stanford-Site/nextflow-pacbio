#!/usr/bin/env nextflow

params.sample = "HG002"
params.threads = "12"
params.mem = "96G"
params.movie_bam = "/oak/stanford/groups/smontgom/jonnguye/old_test_files/HG002/sub_m84011_220902_175841_s1.hifi_reads.bam"
params.reference = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/human_GRCh38_no_alt_analysis_set.fasta"
params.reference_index = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/human_GRCh38_no_alt_analysis_set.fasta.fai"
params.reference_mmi = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/human_GRCh38_no_alt_analysis_set.mmi"
params.reference_name = "GRCh38"
params.outDir = "/oak/stanford/groups/smontgom/jonnguye/test_nextflow/output"
params.reference_tandem_repeat_bed = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/human_GRCh38_no_alt_analysis_set.trf.bed"
params.sex = "MALE"
params.trgt_reference_tandem_repeat_bed = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/trgt/human_GRCh38_no_alt_analysis_set.trgt.v0.3.4.bed"
params.reference_hificnv_exclude_bed = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/hificnv/cnv.excluded_regions.common_50.hg38.bed.gz" 
params.reference_hificnv_exclude_bed_index = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/hificnv/cnv.excluded_regions.common_50.hg38.bed.gz.tbi" 
params.reference_hificnv_expected_bed_female = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/hificnv/female_expected_cn.hg38.bed" 
params.reference_hificnv_expected_bed_male = "/oak/stanford/groups/smontgom/jonnguye/ref/pb-wdl-ref/dataset/GRCh38/hificnv/male_expected_cn.hg38.bed" 

/*
 * Align sample to reference
 */
process pbmm2_align {
    conda 'pbmm2=1.10.0'
    clusterOptions = "--cpus-per-task=12 --mem=32GB --time=1:00:00 --account=smontgom" 
    publishDir params.outDir, mode: 'copy'


    input:
    val sample
    path movie_bam
    val basename
    path reference_mmi
    val reference_name

    output:
    path "${sample}.${basename}.${reference_name}.aligned.bam", emit:alignedBam
    path "${sample}.${basename}.${reference_name}.aligned.bam.bai", emit: alignedBamIndex
 
    """
    pbmm2 align -j 12 -m 32G --sample $sample --log-level INFO --sort --unmapped $reference_mmi $movie_bam ${sample}.${basename}.${reference_name}.aligned.bam
    """
}

process extract_read_length_and_qual {
    conda '/home/jonnguye/micromamba/envs/python3'
    clusterOptions = "--cpus-per-task=1 --mem=8GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'

    input:
    path bam
    val sample_id
    val movie

    output:
    path "${sample_id}.${movie}.read_length_and_quality.tsv"

    """
    /oak/stanford/groups/smontgom/jonnguye/test_nextflow/scripts/extract_read_length_and_qual.py ${bam} > ${sample_id}.${movie}.read_length_and_quality.tsv
    """
}

process pbsv_discover {
    conda 'pbsv=2.9.0'
    clusterOptions = "--cpus-per-task=2 --mem=8GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'

    input:
    path alignedBam
    path alignedBamIndex
    path reference_tandem_repeat_bed

    output:
    path "${alignedBam.getBaseName()}.svsig.gz"

    script:
    """
    pbsv discover --log-level INFO  --hifi --tandem-repeats ${reference_tandem_repeat_bed} ${alignedBam} ${alignedBam.getBaseName()}.svsig.gz
    """
}

process pbsv_call {
    conda 'pbsv=2.9.0'
    clusterOptions="--cpus-per-task=8 --mem=64GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'
    
    input:
    val sampleId
    path svsig
    path reference
    path referenceIndex
    val referenceName

    output:
    path "${sampleId}.${referenceName}.pbsv.vcf", emit: vcf

    """
    pbsv call --hifi --min-sv-length 20 --num-threads ${params.threads} ${reference} ${svsig} ${sampleId}.${referenceName}.pbsv.vcf
    """

}

process zip_index {
    conda "htslib=1.18"
    clusterOptions = "--cpus-per-task=4 --mem=4GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'

    input:
    path vcf

    output:
    path "${vcf.getBaseName()}.vcf.gz", emit: zipped_vcf
    path "${vcf.getBaseName()}.vcf.gz.tbi", emit: zipped_vcf_index

    """
    bgzip --threads 4 --stdout ${vcf} > ${vcf.getBaseName()}.vcf.gz
    tabix --preset vcf ${vcf.getBaseName()}.vcf.gz
    """
}

process deepvariant_make_examples {
    container "/oak/stanford/groups/smontgom/jonnguye/sif/deepvariant.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=1 --mem=8GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    maxForks 8
    
    input:
    //path(input), path(index), path(intervals)
    //#path(fasta)
    //#path(fai)
    //#path(gzi)
    val sample_id
    path aligned_bam
    path bai
    path reference
    path reference_index
    val task_start_index
    val tasks_per_shard
    val total_deepvariant_tasks

    output:
    path "${sample_id}.${task_start_index}.example_tfrecords.tar.gz", emit: example_tfrecords
    path "${sample_id}.${task_start_index}.nonvariant_site_tfrecords.tar.gz", emit: nonvariant_site_tfrecords
    path "example_tfrecords/${sample_id}.examples.tfrecord*${total_deepvariant_tasks}.gz", emit: example_tasks
    path "nonvariant_site_tfrecords/${sample_id}.gvcf.tfrecord*${total_deepvariant_tasks}.gz", emit: nonvariant_tasks

   script:
    """
set -euo pipefail

mkdir example_tfrecords nonvariant_site_tfrecords

echo "DeepVariant version: \$VERSION"

/opt/deepvariant/bin/make_examples \
--norealign_reads \
--vsc_min_fraction_indels 0.12 \
--pileup_image_width 199 \
--track_ref_reads \
--phase_reads \
--partition_size=25000 \
--max_reads_per_partition=600 \
--alt_aligned_pileup=diff_channels \
--add_hp_channel \
--sort_by_haplotypes \
--parse_sam_aux_fields \
--min_mapping_quality=1 \
--mode calling \
--ref ${reference} \
--reads ${aligned_bam} \
--examples example_tfrecords/${sample_id}.examples.tfrecord@${total_deepvariant_tasks}.gz \
--gvcf nonvariant_site_tfrecords/${sample_id}.gvcf.tfrecord@${total_deepvariant_tasks}.gz \
--task ${task_start_index}

echo "Complete"

tar -zcvf ${sample_id}.${task_start_index}.example_tfrecords.tar.gz example_tfrecords
tar -zcvf ${sample_id}.${task_start_index}.nonvariant_site_tfrecords.tar.gz nonvariant_site_tfrecords
"""
}

def get_shard_indices(shards, tasks_per_shard) {
    int counter = 0
    def list = []
    int index = 0
    while (counter < shards){
        
        index = counter * tasks_per_shard;
        list.add(index)
        counter += 1
    }
    return list
}

process deepvariant_call_variants{
    container "/oak/stanford/groups/smontgom/jonnguye/sif/deepvariant.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=16 --mem=64GB --time=4:00:00 --account=smontgom"
    //queue = "batch"
    maxForks 1
    
    input:
    val sample_id
    val reference_name
    path example_tfrecords
    val total_deepvariant_tasks

    output:
    path "${sample_id}.${reference_name}.call_variants_output.tfrecord.gz", emit: tfrecord
   
    script:
	"""
    set -euo pipefail
	deepvariant_model_path="/opt/models/pacbio/model.ckpt"

	echo "DeepVariant version: \$VERSION"

	/opt/deepvariant/bin/call_variants \
		--outfile "${sample_id}.${reference_name}.call_variants_output.tfrecord.gz" \
		--examples "${sample_id}.examples.tfrecord@${total_deepvariant_tasks}.gz" \
		--checkpoint "\${deepvariant_model_path}"
    """
}

process deepvariant_postprocess_variants{
    container "/oak/stanford/groups/smontgom/jonnguye/sif/deepvariant.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=2 --mem=32GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    maxForks 1
    
    input:
    val sample_id
    path tfrecord
    path reference
    path reference_index
    val reference_name
    path nonvariant_tfrecords
    val total_deepvariant_tasks

    output:
    path "${sample_id}.${reference_name}.deepvariant.vcf.gz", emit: vcf
	path "${sample_id}.${reference_name}.deepvariant.vcf.gz.tbi", emit: vcf_index
    path "${sample_id}.${reference_name}.deepvariant.g.vcf.gz", emit: gvcf
	path "${sample_id}.${reference_name}.deepvariant.g.vcf.gz.tbi", emit:gvcf_index

   
    script:
    """
    /opt/deepvariant/bin/postprocess_variants \
			--vcf_stats_report=false \
			--ref ${reference} \
			--infile ${tfrecord} \
			--outfile ${sample_id}.${reference_name}.deepvariant.vcf.gz \
			--nonvariant_site_tfrecord_path "${sample_id}.gvcf.tfrecord@${total_deepvariant_tasks}.gz" \
			--gvcf_outfile ${sample_id}.${reference_name}.deepvariant.g.vcf.gz
	"""

}

process bcftools_on_deepvariant { 
    container "/oak/stanford/groups/smontgom/jonnguye/sif/bcftools.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=2 --mem=8GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    
    input:
    val stats_params
    path reference
    path reference_index
    path vcf
    path vcf_index


    output:
    path "${vcf.simpleName}.vcf.stats.txt", emit: stats
	path "${vcf.simpleName}.bcftools_roh.out", emit: roh_out
	path "${vcf.simpleName}.roh.bed", emit: roh_bed
   
    shell:
    """
set -euo pipefail
bcftools --version

bcftools stats \
--threads 1 \
${stats_params} \
--fasta-ref ${reference} \
${vcf} \
> ${vcf.simpleName}.vcf.stats.txt

echo "STATS"

bcftools roh \
--threads 1 \
--AF-dflt 0.4 \
${vcf} \
> ${vcf.simpleName}.bcftools_roh.out

echo "ROH"

echo -e "#chr\\tstart\\tend\\tqual" > ${vcf.simpleName}.roh.bed
awk -v OFS='\t' '\$1=="RG" {{ print \$3, \$4, \$5, \$8 }}' \
${vcf.simpleName}.bcftools_roh.out \
>> ${vcf.simpleName}.roh.bed

echo "ROH_STATS"

"""

}

process hiphase {
    container "/oak/stanford/groups/smontgom/jonnguye/sif/hiphase.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=16 --mem=64GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    
    input:
    val sample 
    path deepvariant_vcf
    path deepvariant_vcf_index
    path pbsv_vcf
    path pbsv_vcf_index
    path aligned_bam
    path aligned_bam_index
    path reference
    path reference_index
    val reference_name

    output:
    path "${sample}.deepvariant.phased.vcf.gz", emit: deepvariant_phased_vcf
    path "${sample}.deepvariant.phased.vcf.gz.tbi", emit: deepvariant_phased_vcf_index
    path "${sample}.pbsv.phased.vcf.gz", emit: pbsv_phased_vcf 
    path "${sample}.pbsv.phased.vcf.gz.tbi", emit: pbsv_phased_vcf_index 
    path "${sample}.${aligned_bam.baseName}.haplotagged.bam", emit: haplotagged_bam
    path "${sample}.${aligned_bam.baseName}.haplotagged.bam.bai", emit: haplotagged_bam_index
    path "${sample}.${reference_name}.hiphase.stats.tsv" 
    path "${sample}.${reference_name}.hiphase.blocks.tsv"
    
    shell:
    """
    set -euo pipefail

	hiphase --version

    #// for future use with joint calling
    #// String haplotags_param = if length(haplotagged_bam_names) > 0 then "--haplotag-file ~{id}.~{refname}.hiphase.haplotags.tsv" else ""
	
hiphase --threads 16 \
--sample-name ${sample} \
--vcf ${deepvariant_vcf} \
--vcf ${pbsv_vcf} \
--output-vcf "${sample}.deepvariant.phased.vcf.gz" \
--output-vcf "${sample}.pbsv.phased.vcf.gz" \
--bam ${aligned_bam} \
--output-bam "${sample}.${aligned_bam.baseName}.haplotagged.bam" \
--reference ${reference} \
--summary-file ${sample}.${reference_name}.hiphase.stats.tsv \
--blocks-file ${sample}.${reference_name}.hiphase.blocks.tsv \
--global-realignment-cputime 300

bcftools index --tbi --threads 15 "${sample}.deepvariant.phased.vcf.gz" 
bcftools index --tbi --threads 15 "${sample}.pbsv.phased.vcf.gz"

    """
}

process mosdepth {
    container "/oak/stanford/groups/smontgom/jonnguye/sif/mosdepth.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=4 --mem=8GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    
    input:
    path aligned_bam
    path aligned_bam_index

    output:
    path "${aligned_bam.baseName}.mosdepth.summary.txt"
    path "${aligned_bam.baseName}.regions.bed.gz"

    shell:
    """
set -euo pipefail

mosdepth --version

mosdepth \
--threads 3 \
--by 500 \
--no-per-base \
--use-median \
${aligned_bam.baseName} \
${aligned_bam}
    """
}

process trgt {
    container "/oak/stanford/groups/smontgom/jonnguye/sif/trgt.sif"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=4 --mem=8GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    
    input:
    val sex
    path reference
    path reference_index
    path tandem_repeat_bed
    path bam
    path bam_index

    output:
    path "${bam.baseName}.trgt.spanning.bam", emit: spanning_reads 
    path "${bam.baseName}.trgt.vcf.gz", emit: repeat_vcf

    script:
    if (sex == "MALE") {
        karyotype = "XY"
    }
    else {
        karyotype = "XX"
    }
    """
set -euo pipefail

trgt --version

trgt\
 --threads 4 \
 --karyotype ${karyotype} \
 --genome ${reference} \
 --repeats ${tandem_repeat_bed} \
 --reads ${bam} \
 --output-prefix ${bam.baseName}.trgt

"""
}

process sort_trgt_vcf {
    conda "/home/jonnguye/micromamba/envs/common"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=4 --mem=8GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    
    input:
    path vcf 

    output:
    path "*.trgt.sorted.vcf.gz", emit: sorted_vcf 
    path "*.trgt.sorted.vcf.gz.tbi", emit: sorted_vcf_index

    script:
    """
    bcftools --version
    trimmed_name=${vcf}
    trimmed_name=\${trimmed_name%.vcf.gz}
    echo \$trimmed_name
    bcftools sort --output-type z --output "\${trimmed_name}.sorted.vcf.gz" "${vcf}"
    bcftools index --threads 3 --tbi "\${trimmed_name}.sorted.vcf.gz"
"""
}

process sort_trgt_bam {
    conda "/home/jonnguye/micromamba/envs/common"
    publishDir params.outDir, mode: 'copy'
    clusterOptions = "--cpus-per-task=4 --mem=8GB --time=2:00:00 --account=smontgom"
    //queue = "batch"
    
    input:
    path bam 

    output:
    path "${bam.baseName}.sorted.bam", emit: sorted_bam
    path "${bam.baseName}.sorted.bam.bai", emit: sorted_bam_index

    script:
    """
samtools --version
samtools sort -@ 3 -o "${bam.baseName}.sorted.bam" ${bam}
samtools index -@ 3 "${bam.baseName}.sorted.bam"
    """
}


process coverage_dropouts {
    conda '/home/jonnguye/micromamba/envs/python3'
    clusterOptions = "--cpus-per-task=2 --mem=8GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'

    input:
    path bam
    path bam_index
    path trgt_reference_tandem_repeat_bed
    val sample
    val reference_name

    output:
    path "${sample}.${reference_name}.trgt.dropouts.txt", emit: coverage_dropouts
    
    script:
    """
set -euo pipefail
/oak/stanford/groups/smontgom/jonnguye/test_nextflow/scripts/check_trgt_coverage.py \
    ${trgt_reference_tandem_repeat_bed} \
	${bam} \
	> "${sample}.${reference_name}.trgt.dropouts.txt"
    """
}

process cpg_pileup {
    container "/oak/stanford/groups/smontgom/jonnguye/sif/pb-cpg-tools.sif"
    clusterOptions = "--cpus-per-task=12 --mem=48GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'
    
    input:
    path bam
    path bam_index
    val sample
    val reference_name
    path reference
    path reference_index

    output:
    path "${sample}.${reference_name}.*.bed", emit: pileup_beds
    path "${sample}.${reference_name}.*.bw", emit: pileup_bigwigs

    script:
    """
	set -euo pipefail

	aligned_bam_to_cpg_scores --version

	aligned_bam_to_cpg_scores \
    --threads 12 \
	--bam ${bam} \
	--ref ${reference} \
	--output-prefix "${sample}.${reference_name}"  \
	--min-mapq 1 \
	--min-coverage 10 \
    --model "\${PILEUP_MODEL_DIR}/pileup_calling_model.v1.tflite"
   """ 
}

process paraphase {
    container "/oak/stanford/groups/smontgom/jonnguye/sif/paraphase.sif"
    clusterOptions = "--cpus-per-task=4 --mem=8GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'

    input:
    val sample
    path bam
    path bam_index
    path reference
    path reference_index

    output:
    path "${sample}.paraphase/${sample}.json"
    path "${sample}.paraphase/${sample}_realigned_tagged.bam"
    path "${sample}.paraphase/${sample}_realigned_tagged.bam.bai"
    path "${sample}.paraphase/${sample}_vcfs/*.vcf", optional: true

    script:
	"""
    set -euo pipefail

	paraphase --version

	paraphase \
		--threads 4 \
		--bam ${bam} \
		--reference ${reference} \
		--out "${sample}.paraphase"
    """
}

process hificnv{
    container "/oak/stanford/groups/smontgom/jonnguye/sif/hificnv.sif"
    clusterOptions = "--cpus-per-task=8 --mem=16GB --time=1:00:00 --account=smontgom"
    publishDir params.outDir, mode: 'copy'

    input:
    val sample
    val sex
    path bam
    path bam_index
    path phased_vcf
    path phased_vcf_index
    path reference
    path reference_index
    path exclude_bed
    path exclude_bed_index
    path expected_bed_male
    path expected_bed_female

    output:
    path "hificnv.${sample}.vcf.gz"
    path "hificnv.${sample}.vcf.gz.tbi"
    path "hificnv.${sample}.copynum.bedgraph"
    path "hificnv.${sample}.depth.bw"
    path "hificnv.${sample}.maf.bw"

    script:
    if (sex == "MALE") {
        expected_bed = "${expected_bed_male}"
    }
    else {
        expected_bed = "${expected_bed_female}"
    }

    """
    set -euo pipefail

	hificnv --version

	hificnv \
		--threads 8 \
		--bam ${bam} \
		--ref ${reference} \
		--maf ${phased_vcf} \
		--exclude ${exclude_bed} \
		--expected-cn ${expected_bed} \
		--output-prefix hificnv 

	bcftools index --tbi "hificnv.${sample}.vcf.gz"
    """
}


/*
 * Define the workflow
 */
workflow {
    basename = file(params.movie_bam).getBaseName()
    alignedBam_ch = pbmm2_align(params.sample, params.movie_bam, basename, params.reference_mmi, params.reference_name) 
    extract_read_length_and_qual(params.movie_bam, params.sample, basename)
    svsig_ch = pbsv_discover(alignedBam_ch.alignedBam, alignedBam_ch.alignedBamIndex, params.reference_tandem_repeat_bed)
    pbsv_vcf_ch = pbsv_call(params.sample, svsig_ch, params.reference, params.reference_index, params.reference_name)
    pbsv_vcf_zip_ch = zip_index(pbsv_vcf_ch)
    deepvariant_tasks = 64
    tasks_per_shard = 1 
    shard_indices = Channel.from(get_shard_indices(deepvariant_tasks, tasks_per_shard))
    deepvariant_make_examples_ch = deepvariant_make_examples(params.sample, alignedBam_ch.alignedBam, alignedBam_ch.alignedBamIndex, params.reference, params.reference_index, shard_indices, tasks_per_shard, deepvariant_tasks)
    example_tasks = deepvariant_make_examples_ch.example_tasks.collect()
    nonvariant_tasks = deepvariant_make_examples_ch.nonvariant_tasks.collect()
    deepvariant_calls_ch = deepvariant_call_variants(params.sample, params.reference_name, example_tasks, deepvariant_tasks)
    deepvariant_postprocess_variants_ch = deepvariant_postprocess_variants(params.sample, deepvariant_calls_ch.tfrecord, params.reference, params.reference_index, params.reference_name, nonvariant_tasks, deepvariant_tasks)
    
    stats_params = "--apply-filters PASS --samples ${params.sample}"
    bcftools_ch = bcftools_on_deepvariant(stats_params, params.reference, params.reference_index, deepvariant_postprocess_variants_ch.vcf, deepvariant_postprocess_variants_ch.vcf_index)

    hiphase_ch = hiphase(params.sample, 
        deepvariant_postprocess_variants_ch.vcf, deepvariant_postprocess_variants_ch.vcf_index,
        pbsv_vcf_zip_ch.zipped_vcf , pbsv_vcf_zip_ch.zipped_vcf_index,
        alignedBam_ch.alignedBam, alignedBam_ch.alignedBamIndex,
        params.reference, params.reference_index, params.reference_name
    )

    mosdepth_ch = mosdepth(hiphase_ch.haplotagged_bam, hiphase_ch.haplotagged_bam_index)
    trgt_ch = trgt(params.sex, params.reference, params.reference_index, params.trgt_reference_tandem_repeat_bed, hiphase_ch.haplotagged_bam, hiphase_ch.haplotagged_bam_index)
    trgt_sort_vcf_ch = sort_trgt_vcf(trgt_ch.repeat_vcf)
    trgt_sort_bam_ch = sort_trgt_bam(trgt_ch.spanning_reads)
    coverage_dropouts(hiphase_ch.haplotagged_bam,hiphase_ch.haplotagged_bam_index, params.trgt_reference_tandem_repeat_bed,
        params.sample, params.reference_name)
    cpg_pileup(hiphase_ch.haplotagged_bam, hiphase_ch.haplotagged_bam_index, params.sample, params.reference_name, 
        params.reference, params.reference_index)    
    paraphase(params.sample, hiphase_ch.haplotagged_bam, hiphase_ch.haplotagged_bam_index,
        params.reference, params.reference_index)
    hificnv(params.sample, params.sex, hiphase_ch.haplotagged_bam, hiphase_ch.haplotagged_bam_index,
        hiphase_ch.deepvariant_phased_vcf, hiphase_ch.deepvariant_phased_vcf_index,
        params.reference, params.reference_index,
        params.reference_hificnv_exclude_bed, params.reference_hificnv_exclude_bed_index,
        params.reference_hificnv_expected_bed_male, params.reference_hificnv_expected_bed_female
    )
}