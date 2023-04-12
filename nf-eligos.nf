#!/usr/bin/env nextflow
/*
========================================================================================
                         MaestSi/nf-eligos
========================================================================================
 MaestSi/nf-eligos analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/nf-eligos
----------------------------------------------------------------------------------------
*/

def helpMessage() {
		log.info"""
	Usage:
	nextflow -c nf-eligos.conf run nf-eligos.nf --samples="/path/to/samples.txt" --resultsDir="/path/to/resultsDir" -profile docker

	Mandatory argument:
	-profile                                                 Configuration profile to use. Available: docker, singularity
	Other mandatory arguments which may be specified in the nf-eligos.conf file

	--samples                                                Path to the tab-separated sample file including sample name, condition and path to fastq file
        --baseline_condition                                     Condition to be considered as the baseline, must match one of the conditions in the samples file
        --min_depth                                              Minimum number of reads
        --max_depth                                              Maximum number of reads
        --spliced_alignment_flag                                 Flag for splice-aware alignment, set to true for genome alignment and to false for transcriptome alignment
	--min_mapq                                               Minimum mapping quality for filtering alignments
        --resultsDir                                             Path to a folder where to store results
	--reference_fasta                                        Path to the reference fasta file
	--bed_file                                               Path to regions of interest bed file
        --pval_thr                                               p-value threshold
        --adjPval_thr                                            adjusted p-value threshold
        --oddR_thr                                               Odds Ratio threshold
        --esb_thr                                                Threshold on %Error of Specific Bases to be considered for de novo motifs discovery
        --sb                                                     Selected basis for filtering modification of interest
        --opt_args                                               Other optional arguments (e.g. "-bcf <file.bcf> -m <model.json>")
	""".stripIndent()
}

// Show help message
if (params.help) {
	helpMessage()
	exit 0
}

// Input of reference fasta
Channel
	.fromPath(params.reference_fasta, checkIfExists:true)
	.into{reference_fasta_minimap2; reference_fasta_eligosPair; reference_fasta_eligosRbem}

// Input of bed file
Channel
        .fromPath(params.bed_file, checkIfExists:true)
        .into{bed_file_eligosPair; bed_file_eligosRbem}

// Input of sample names, conditions and FASTQ path
Channel
	.fromPath( params.samples )
	.splitCsv(header: true, sep:'\t')
	.map{ row-> tuple(row.sample, row.condition, file(row.fastq)) }
	.set{ samples_minimap2 }

// Reference alignment
process minimap2 {
	input:
	tuple val(sample), val(condition), val(fastq) from samples_minimap2

	each file('reference.fasta') from reference_fasta_minimap2
		
	output:
	tuple val(condition), val(sample), file('minimap.filtered.bam') into minimap2_bamMerge

    script:
    if(params.minimap2)
    """
        mkdir -p ${params.resultsDir}/${condition}/${sample}/Alignment

        if ${params.spliced_alignment_flag} ; then 
		minimap2 -ax splice -k14 -t ${task.cpus} reference.fasta ${fastq} | samtools view -hSb | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.bam
		samtools view ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.bam -bh -q ${params.min_mapq} -F 2308 | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam
	else
		minimap2 -ax map-ont -k14 -t ${task.cpus} reference.fasta ${fastq} | samtools view -hSb | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.bam
                samtools view ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.bam -bh -q ${params.min_mapq} -F 2324 | samtools sort -@ ${task.cpus} -o ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam
		samtools index -@ ${task.cpus} ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam
	fi

	ln -s ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam ./minimap.filtered.bam
	ln -s ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam.bai ./minimap.filtered.bam.bai
	"""
	else
	"""
	ln -s ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam ./minimap.filtered.bam
        ln -s ${params.resultsDir}/${condition}/${sample}/Alignment/minimap.filtered.bam.bai ./minimap.filtered.bam.bai		
	"""
}

// From a single channel for all the alignments to one channel for each condition
minimap2_bamMerge.groupTuple(by:0)
.set {minimap2_bamMerge_grouped}

process bamMerge {
	input:
	tuple val(condition), val(sample), file('minimap.filtered.*.bam') from minimap2_bamMerge_grouped

	output:
	tuple val(condition), val(sample) into bamMerge_eligosPairTmp
	tuple val(condition), val(sample) into bamMerge_eligosRbemTmp

    script:
    if(params.bamMerge)
	"""
		mkdir -p ${params.resultsDir}/${condition}/Alignment/
		samtools merge -f ${params.resultsDir}/${condition}/Alignment/${condition}.filt.bam minimap.filtered.*.bam
		samtools index ${params.resultsDir}/${condition}/Alignment/${condition}.filt.bam
	"""
	else
	"""
		ln -s ${params.resultsDir}/${condition}/Alignment/${condition}.filt.bam ./${condition}.filt.bam
		ln -s ${params.resultsDir}/${condition}/Alignment/${condition}.filt.bam.bai ./${condition}.filt.bam.bai
	"""
}

// From a single channel for all the alignments to one channel for each condition
bamMerge_eligosPairBaseline=Channel.create()
bamMerge_eligosPairOther=Channel.create()

bamMerge_eligosPairTmp
.choice( bamMerge_eligosPairBaseline, bamMerge_eligosPairOther ) { a -> a[0] == params.baseline_condition ? 0 : 1 } 

process eligosPair {
        input:
        tuple val('conditionBaseline'), val('sample') from bamMerge_eligosPairBaseline.collect()
        tuple val('conditionTest'), val('sample') from bamMerge_eligosPairOther
        each file('file.bed') from bed_file_eligosPair
        each file('reference.fasta') from reference_fasta_eligosPair
        output:

    script:
    if(params.eligosPair)
        """
                mkdir -p ${params.resultsDir}/${conditionTest}/eligosPair/
                
                eligos2 pair_diff_mod -t ${task.cpus} \
		-tbam ${params.resultsDir}/${conditionTest}/Alignment/${conditionTest}.filt.bam \
		-cbam ${params.resultsDir}/${conditionBaseline}/Alignment/${conditionBaseline}.filt.bam \
		-reg file.bed \
		-ref reference.fasta \
                --min_depth ${params.min_depth} --max_depth ${params.max_depth} \
                --pval ${params.pval_thr} --oddR ${params.oddR_thr} --esb ${params.esb_thr} --adjPval ${params.adjPval_thr} \
		-o ${params.resultsDir}/${conditionTest}/eligosPair ${params.opt_args}

                res=\$(find ${params.resultsDir}/${conditionTest}/eligosPair | grep "combine\\.txt")
                eligos2 filter -i \$res -sb ${params.sb} --homopolymer --oddR ${params.oddR_thr} --esb ${params.esb_thr} --adjPval ${params.adjPval_thr}
                
                res_filt=\$(find ${params.resultsDir}/${conditionTest}/eligosPair | grep "filtered\\.txt")
                eligos2 bedgraph -i \$res_filt -sb ${params.sb} --signal oddR --homopolymer
        """
	else
	"""
		echo "Skipped"
	"""
}

// From a single channel for all the alignments to one channel for each condition
bamMerge_eligosRbemTmp.groupTuple(by:0)
.set { bamMerge_eligosRbem }

process eligosRbem {
        input:
        tuple val(condition), val(sample) from bamMerge_eligosRbem
        each file('file.bed') from bed_file_eligosRbem
        each file('reference.fasta') from reference_fasta_eligosRbem
        output:
                
    script:
    if(params.eligosRbem)
        """
                mkdir -p ${params.resultsDir}/${condition}/eligosRbem/
                
                eligos2 rna_mod -t ${task.cpus} \
                -i ${params.resultsDir}/${condition}/Alignment/${condition}.filt.bam \
                -reg file.bed \
                -ref reference.fasta \
                --pval ${params.pval_thr} --oddR ${params.oddR_thr} --esb ${params.esb_thr} --adjPval ${params.adjPval_thr} \
                --min_depth ${params.min_depth} --max_depth ${params.max_depth} \
                -o ${params.resultsDir}/${condition}/eligosRbem ${params.opt_args}
               
                res=\$(find ${params.resultsDir}/${condition}/eligosRbem | grep "Ext0\\.txt")
                eligos2 bedgraph -i \$res -sb ${params.sb} --signal ESB --homopolymer
        """
	else
	"""
		echo "Skipped"
	"""
}
