# star-rsem-scrnaseq-nf
`star-rsem-scRNAseq`: a workflow for reprocessing of Maynard et al. (Cell 2020, https://doi.org/10.1016/j.cell.2020.07.017)  
`STAR` parameters are taken from BAM header of author provided BAM file  

# Usage
Point `star-rsem-scrnaseq-nf` at a directory containing fastqs split by cell barcode to generate alignments and RSEM counts  

```
nextflow run star-rsem-scrnaseq-nf --reads [full path] --genome [full path] --rsem_index [full path] --outdir [full path]
```

By default `star-rsem-scrnaseq-nf` will look for `--genome`, and `--rsem_index` in HCI Core group space on Redwood

## Required
+ `--reads`
    + Full path to directory containing fastqs demultiplexed by cell barcode. Assumes fastqs are gzipped. Default is `./fastq`
+ `--genome`
    + Full path to STAR formatted reference. Defaults to Hg38 in HCI Core group space.  
+ `--rsem_index`
    + Full path to RSEM index. Defaults to Hg38 in HCI Core group space
+ `--outdir`
    + Path to publish results.  

## Execution
Standard profile submits jobs to `slurm` scheduler via the shared queue.    

## Output
`star-rsem-scrnaseq-nf` will publish BAMs and indicies to `--outdir/bams` and RSEM results to `--outdir/rsem_out`.  
Nextflow timeline, report, and trace are published to `./logs`

## Versions (by CHPC modules)
`STAR 2.7.8a`  
`samtools 1.12`  
`rsem 1.3.0`  
