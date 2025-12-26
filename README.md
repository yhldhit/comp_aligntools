# Running the experiment requires the following steps:

1. Prepare the data, including raw sequencing reads (FASTQ.gz files) and the reference genome (e.g., hg38 for human samples).
2. Set up the environment by installing toolkits such as `cutadapt`, `Bismark`, `samtools`, and `Picard`; modify `bismark_pipeline.sh` and run the Bismark alignment.  
3. Run Parabricks following the reference: [https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam_meth.html](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam_meth.html); prepare the reference genome in `bwa-meth` format, modify and run `parabricks_pipeline.sh` to perform Parabricks alignment.  
4. Use any tool (e.g., `wgbstools`) to extract corresponding BED files from the generated BAM files.  
5. Use the script `compbed_venn.py` to compare two BED files.
