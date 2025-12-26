#!/bin/bash
#SBATCH --mail-user=1317575996@qq.com
#SBATCH --mail-type=ALL
#SBATCH -N 3
#SBATCH --mem=300GB
#SBATCH --ntasks-per-node=50
#SBATCH --output=Job_bismark.%j.out
#SBATCH --error=Job_bismark.%j.err
#SBATCH --partition=cpu
source ~/anaconda3/bin/activate
conda activate py38
python /home/lijia/yanhongliang/proj/comp_aligntools/src/bismark_pipeline.py -i ./softlink_fastq -o ./softlink_fastq_result -p 15 --bismark_reference ./hg38_gene --reference_fasta /home/wang_yanni/genome_file/hg38.fa --chrom_size_path /home/wang_yanni/genome_file/hg38.chrom.sizes
