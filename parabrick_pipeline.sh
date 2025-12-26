#!/bin/bash
#SBATCH --job-name=fq2bam_meth_cutadapt
#SBATCH --output=logs/fq2bam_meth_%j.out
#SBATCH --error=logs/fq2bam_meth_%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:2
#SBATCH --nodelist=gpu01
#SBATCH --cpus-per-task=60
#SBATCH --mem=100G
#SBATCH --time=96:00:00

# =====================================================
# 环境与初始化
# =====================================================
set -e
mkdir -p logs out trimmed tmp

# 镜像与参考文件路径
IMAGE_NAME="nvcr.io/nvidia/clara/clara-parabricks:4.5.1-1"
IMAGE_TAR="/home/yan_hongliang/proj/parabrick2bam/parabricks_4.5.1.tar"
FASTQ_DIR="/home/yan_hongliang/proj/parabrick2bam/demo1225/fastq"
REFERENCE_FILE="hg38_bwa_meth/hg38.fa"

# Cutadapt 参数（与 Bismark 一致）
CPU=60
MIN_LENGTH=30
QUALITY=20
OVERLAP=6
R1_ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
R2_ADAPTER="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# =====================================================
# Step 1. 加载镜像
# =====================================================
echo "[$(date)] Checking Parabricks image..."
if ! docker images | grep -q "clara-parabricks"; then
    if [ -f "$IMAGE_TAR" ]; then
        echo "Image not found locally, loading from TAR..."
        docker load -i "$IMAGE_TAR"
    else
        echo "ERROR: Cannot find $IMAGE_TAR"
        exit 1
    fi
fi

# 激活 Conda 环境（Cutadapt 与 Picard）
source /home/wang_yanni/anaconda3/bin/activate
conda activate py38

# =====================================================
# Step 2. 遍历 FASTQ 样本
# =====================================================
for fq1 in ${FASTQ_DIR}/*.R1.fastq.gz; do
    CELL_ID=$(basename "$fq1" .R1.fastq.gz)
    fq2="${FASTQ_DIR}/${CELL_ID}.R2.fastq.gz"

    if [ ! -f "$fq2" ]; then
        echo "⚠️ Missing R2 for $CELL_ID, skipping..."
        continue
    fi

    DEDUP_BAM="out/${CELL_ID}.dedup.bam"
    if [ -f "$DEDUP_BAM" ]; then
        echo "✅ Skipping $CELL_ID (dedup.bam already exists)"
        continue
    fi

    echo "[$(date)] Processing $CELL_ID ..."

    TRIM_R1="trimmed/${CELL_ID}_R1.trimmed.fq.gz"
    TRIM_R2="trimmed/${CELL_ID}_R2.trimmed.fq.gz"
    FINAL_BAM="out/${CELL_ID}.final.bam"
    FILTERED_BAM="out/${CELL_ID}.final.filtered.bam"
    TMP_DIR="tmp/${CELL_ID}"
    mkdir -p "$TMP_DIR"

    # =================================================
    # Step 3. Cutadapt trimming
    # =================================================
    echo "[$(date)] Running Cutadapt trimming..."
    cutadapt -j ${CPU} -a ${R1_ADAPTER} -A ${R2_ADAPTER} \
        -q ${QUALITY} -m ${MIN_LENGTH} -O ${OVERLAP} \
        -o ${TRIM_R1} -p ${TRIM_R2} \
        ${fq1} ${fq2} > trimmed/${CELL_ID}_cutadapt.log 2>&1

    # =================================================
    # Step 4. Parabricks 比对
    # =================================================
    echo "[$(date)] Running Parabricks fq2bam_meth..."
    docker run --rm --gpus '"device=6,7"' \
        -v $(pwd):/workdir \
        --workdir /workdir \
        "$IMAGE_NAME" \
        pbrun fq2bam_meth \
            --ref /workdir/${REFERENCE_FILE} \
            --in-fq /workdir/${TRIM_R1} /workdir/${TRIM_R2} \
            --out-bam /workdir/${FINAL_BAM} \
            --filter-flag 3852 \
            --min-read-length ${MIN_LENGTH} \
            --max-read-length 480 \
            --num-gpus 2 \
            --bwa-nstreams 2 \
            --bwa-cpu-thread-pool 16 \
            --low-memory \
            --tmp-dir /workdir/${TMP_DIR}

    # =================================================
    # Step 5. MAPQ>=10 + QC pass
    # =================================================
    echo "[$(date)] Filtering QC failed + low MAPQ reads..."
    samtools view -@ 8 -b -q 10 -F 512 -o ${FILTERED_BAM} ${FINAL_BAM}
    samtools index ${FILTERED_BAM}

    # =================================================
    # Step 6. 去重
    # =================================================
    echo "[$(date)] Removing PCR duplicates..."
    java -Xmx30g -jar /home/wang_yanni/anaconda3/envs/py38/share/picard-3.4.0-0/picard.jar MarkDuplicates \
        I=${FILTERED_BAM} \
        O=${DEDUP_BAM} \
        M=out/${CELL_ID}.dedup.metrics.txt \
        REMOVE_DUPLICATES=true \
        TMP_DIR=${TMP_DIR}
    samtools index ${DEDUP_BAM}

    # =================================================
    # Step 7. 清理中间文件
    # =================================================
    echo "[$(date)] Cleaning intermediate files for ${CELL_ID}..."
    rm -f ${FINAL_BAM} ${FILTERED_BAM} ${TRIM_R1} ${TRIM_R2}
    rm -rf ${TMP_DIR}

    echo "[$(date)] Finished ${CELL_ID} ✅"
done

echo "[$(date)] All samples processed successfully!"
