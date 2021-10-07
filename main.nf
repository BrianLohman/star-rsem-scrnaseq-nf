#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
if (params.help) {
    log.info"""
    -----------------------------------------------------------------------------
    star-rsem-scRNAseq: a workflow for reprocessing of Maynard et al. (Cell 2020)
    =============================================================================

    Required arguments:
    -------------------

    --reads         Full path to directory with reads. Default: ./fastq

    --outdir        Path to publish results. Default: ./results

    --genome        Full path to STAR reference. Default: Hg38 in core group space

    --rsem_index    Full path to RSEM index. Default: Hg38 in core group space

    
    Description:
    ------------
    Runs STAR with parameters from Maynard et al 2020 (Cell: https://doi.org/10.1016/j.cell.2020.07.017)
    Extracted STAR parameters from BAM file provided by Wei Wu via S3

    Index BAMS (BAMS and indicies are published to params.outdir/bams)

    Run RSEM to get gene and isoform counts (saved to params.outdir/rsem_out)
        NB: Deviates from original publication which use htseq-counts

    -----------------------------------------------------------------------------
    """.stripIndent()
    exit 0
}

// Set default required params
params.reads = "./fastq/*_{1,2}.fastq.gz"
params.genome =  '/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release104/star125'
params.rsem_index = "/uufs/chpc.utah.edu/common/HIPAA/hci-bioinformatics1/atlatl/data/Human/GRCh38/release104/rsem/RSEM"
params.outdir = './results'

// Logging
log.info("\n")
log.info("Reads directory (--reads)      :${params.reads}")
log.info("Genome          (--genome)     :${params.genome}")
log.info("RSEM index      (--rsem_index) :${params.rsem_index}")
log.info("Results         (--outdir)     :${params.outdir}")

// Read pair channel: tuple of pair_id and paired fastqs
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

// Maps each read-pair by using STAR
process star {
    module 'star/2.7.8a'
    tag "$pair_id"
      
    input:
      path(genome)
      tuple val(pair_id), path(reads)
  
    output:
      tuple val(pair_id), path("${pair_id}.bam"), emit: bams
      tuple val(pair_id), path("${pair_id}_aligned_to_transcriptome.out.bam"), emit: rsem_input

    script:
      """
      STAR --genomeDir $genome \
      --genomeLoad LoadAndKeep \
      --readFilesIn $reads \
      --readFilesCommand zcat \
      --runThreadN 12 \
      --outSAMtype BAM SortedByCoordinate \
      --limitBAMsortRAM 34359720776 \
      --outBAMsortingBinsN 100 \
      --quantMode TranscriptomeSAM \
      --outFilterType BySJout \
      --outFilterMultimapNmax 20 \
      --alignSJoverhangMin 8 \
      --alignSJDBoverhangMin 1 \
      --outFilterMismatchNmax 999 \
      --outFilterMismatchNoverLmax 0.04 \
      --alignIntronMin 20 \
      --alignIntronMax 1000000 \
      --alignMatesGapMax 1000000 \
      --outSAMstrandField intronMotif \
      --outSAMattributes NH HI NM MD
      mv Aligned.sortedByCoord.out.bam ${pair_id}.bam
      mv Aligned.toTranscriptome.out.bam ${pair_id}_aligned_to_transcriptome.out.bam
    """
}

// Index BAMs
process index {
    module 'samtools/1.12'
    tag "${pair_id}"
    publishDir "${params.outdir}/bams", mode:"copy"

    input:
      tuple val(pair_id), path(bam)

    output:
      path("${pair_id}.bam")
      path("${pair_id}.bam.bai")    

    script:
      """
      samtools index ${pair_id}.bam
      """
}

// RSEM to get gene and isoform counts
process rsem {
    module 'rsem/1.3.0'
    tag "$pair_id"
    publishDir "${params.outdir}/rsem_out", mode:"copy"
 
    input:
      path(params.rsem_index)
      tuple val(pair_id), path(rsem_bam)
 
    output:
      path("${pair_id}.genes.results")
      path("${pair_id}.isoforms.results")
 
    script:
      """
      rsem-calculate-expression --paired-end -num-threads 8 --alignments --strandedness none --no-bam-output \
      ${rsem_bam} ${params.rsem_index} ${pair_id}
      """
}

workflow {
    star(params.genome, read_pairs_ch)
    index(star.out.bams)
    rsem(params.rsem_index, star.out.rsem_input)
}
