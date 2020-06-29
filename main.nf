//  example:
// nextflow run main.nf --reads "/run/media/alex/Seagate Backup Plus Drive/cincinatti_exome_data/FastQ/09KTF_TTAGGC_L002_R{1,2}_001.fastq.gz" --samplename 09KTF_TTAGGC
 
 // Define the default parameters
params.reads = ""
params.outdir = "./results"
params.samplename = ""
params.readgroup = "${params.samplename}"
params.reference = "~/hg38/Homo_sapiens_assembly38.fasta"

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 


process map_and_sort {
	publishDir "${params.outdir}/${params.samplename}"

	input:
	tuple val(pair_id), path(reads) from read_pairs_ch

	output:
	file "${params.samplename}.sorted.bam" into bamfile_ch
	file "${params.samplename}.sorted.bai" into bamindex_ch
	
	"""
	bwa mem -t 6 -M -R '@RG\\tID:${params.readgroup}\\tSM:${params.samplename}\\tPL:ILLUMINA' $params.reference $reads \
	| gatk SortSam -I /dev/stdin -O ${params.samplename}.sorted.bam --SORT_ORDER=coordinate --CREATE_INDEX=true
	"""	
}

process mark_duplicates {
	publishDir "${params.outdir}/${params.samplename}"

	input:
	file bamfile from bamfile_ch

	output:
	file "${params.samplename}.sorted.dedup.bam" into dedup_ch1, dedup_ch2
	file "${params.samplename}.sorted.dedup.metrics.txt" into dudup_metrics_ch

	script:
	"""
	gatk MarkDuplicates \
	-I $bamfile \
	-O ${params.samplename}.sorted.dedup.bam \
	-M ${params.samplename}.sorted.dedup.metrics.txt
	"""
}


process generate_bqsr {
	publishDir "${params.outdir}/${params.samplename}"

	input:
	file dedup_file from dedup_ch1

	output:
	file "recal_${params.samplename}.txt" into recal_ch

	script:
	"""
	gatk BaseRecalibrator --input $dedup_file --output recal_${params.samplename}.txt \
	--known-sites ~/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
	--reference ${params.reference}
	"""
}

process apply_bqsr {
	publishDir "${params.outdir}/${params.samplename}"

	input:
	file dedup_file from dedup_ch2
	file recal_file from recal_ch

	output:
	file "${params.samplename}.sorted.dedup.recal.bam" into dedup_recal_ch

	script:
	"""
	gatk ApplyBQSR -R ${params.reference} \
	-I=$dedup_file --bqsr-recal-file=$recal_file \
	-O=${params.samplename}.sorted.dedup.recal.bam
	"""
}

process call_haplotypes {
	publishDir "${params.outdir}/${params.samplename}"

	input:
	file dedup_recal_file from dedup_recal_ch

	output:
	file "result_${params.samplename}.vcf" into vcf_ch
	file "${params.samplename}_bamout.bam" into bamout_ch

	script:
	"""
	gatk HaplotypeCaller --input $dedup_recal_file --output result_${params.samplename}.vcf \
	--reference ${params.reference} --bamout ${params.samplename}_bamout.bam
	"""
}

process cnn_score_variants {
	publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

	input:
	file vcf from vcf_ch
	file bamout from bamout_ch

	output:
	file "result_cnn2d_${params.samplename}.vcf" into vcf_cnn2d_ch


	script:
	"""
	gatk CNNScoreVariants \
		-R ${params.reference} \
		-I $bamout \
		-V $vcf \
		-O result_cnn2d_${params.samplename}.vcf \
		--tensor-type read_tensor
	"""
}


process filter_variants {
	publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

	input:
	file vcf_cnn2d from vcf_cnn2d_ch
	

	output:
	file "result_cnn2d_filtered_${params.samplename}.vcf" into vcf_cnn2d_filtered_ch


	script:
	"""
	gatk FilterVariantTranches \
	  -V $vcf_cnn2d \
	  -O result_cnn2d_filtered_${params.samplename}.vcf \
	  --resource ~/hg38/1000G_omni2.5.hg38.vcf.gz \
	  --resource ~/hg38/hapmap_3.3.hg38.vcf.gz \
	  --info-key CNN_2D \
	  --snp-tranche 99.9 \
	  --indel-tranche 95.0
	"""
}

process annotate_csq {
	publishDir "${params.outdir}/${params.samplename}", mode: 'copy'

	input:
	file vcf_cnn2d_filtered from vcf_cnn2d_filtered_ch
	

	output:
	file "result_cnn2d_filtered_csq_${params.samplename}.bcf" into vcf_cnn2d_filtered_csq_ch


	script:
	"""
	bcftools csq -f ${params.reference} \
	 -g ~/hg38/Homo_sapiens.GRCh38.100.gff3.gz \
	  $vcf_cnn2d_filtered -Ob -o result_cnn2d_filtered_csq_${params.samplename}.bcf
	"""
}

