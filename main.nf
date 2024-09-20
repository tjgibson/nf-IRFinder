#! /usr/bin/env nextflow

log.info """
	IRFinder pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()

/* input files:	
 * samplesheet
 * merged, sorted and indexed bam files for each RNA-seq sample
 * fasta file
 * gtf file
 */
 
 /*
 * Build IRFinder reference
 */


 process build_ref {
	publishDir "${params.results_dir}/IRFinder_ref", mode: 'copy'
	
	input:
	path(fasta)
	path(gtf)
	
	output:
	path("IRFinder_ref/${params.genome_name}")
	
	script:
    """
	mkdir -p IRFinder_ref/${params.genome_name}
	
	gunzip -c $fasta > IRFinder_ref/${params.genome_name}/genome.fa
	gunzip -c $gtf > IRFinder_ref/${params.genome_name}/transcripts.gtf
	
	
	IRFinder BuildRefProcess -r IRFinder_ref/${params.genome_name}
	"""
    
    stub:
    """
    mkdir -p IRFinder_ref/${params.genome_name}
    touch  IRFinder_ref/${params.genome_name}/genome.fa
    touch IRFinder_ref/${params.genome_name}/transcripts.gtf
    """	
}

/*
 * Run IRFinder on BAM files
 */

process quantify_IR {
	tag "$meta.sample"
	publishDir "${params.results_dir}/IR_quantification/", mode: 'copy'
	
	input:
	tuple val(meta), path(bam), path(bam_index)
	path(IRFinder_ref)

	
	output:
	tuple val(meta), path("${meta.sample}")
	
	script:
    """
	IRFinder BAM \
        -r $IRFinder_ref \
        -d ${meta.sample} \
        $bam
	"""
	
	stub:
	"""
	mkdir ${meta.sample}
	touch ${meta.sample}/IRFinder-IR-dir.txt
	touch ${meta.sample}/IRFinder-IR-nondir.txt
	
	"""
}

/*
 * Run workflow
 */
 
workflow {
	
	bam_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
	| splitCsv( header:true )
    | map { row ->
        bam_meta = row.subMap('sample')
        [
        	bam_meta, 
        	file(row.reads, checkIfExists: true),
            file(row.read_index, checkIfExists: true)]
    }
	
	ref_ch = build_ref(
	file("${launchDir}/${params.fasta}", checkIfExists: true),
	file("${launchDir}/${params.gtf}", checkIfExists: true)
	)
	
	quantify_IR(
	bam_ch,
	ref_ch
	)	
}