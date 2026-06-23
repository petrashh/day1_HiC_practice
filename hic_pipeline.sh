#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

THREADS=${THREADS:-4}
BASE_URL="https://genedev.bionet.nsc.ru/ftp/_RawReads/2025-05-23MyGenetics"
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
MIN_LEN=70
QUAL_CUTOFF=20
ENZYME="DpnII"
GENOME_SHORT="T2T_human"
GENOME_FA="data/reference/T2T_human.fna"
CHROM_SIZES="data/reference/chrom.sizes"
SITE_POSITIONS="data/reference/restriction_sites_DpnII.txt"
JUICER_DIR="$(pwd)/tools/juicer"

declare -A SAMPLES
SAMPLES=(
  ["MoPh11"]="Copy%20of%20MoPh11_S86_L001"
  ["MoPh14"]="Copy%20of%20MoPh14_S87_L001"
  ["MoPh15"]="Copy%20of%20MoPh15_S88_L001"
)

LOG_DIR="logs"
mkdir -p "${LOG_DIR}"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/pipeline.log"
}

log "Checking prerequisites"
for tool in fastqc cutadapt bwa samtools java wget; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Missing $tool"; exit 1; }
done

for f in "${GENOME_FA}" "${GENOME_FA}.bwt" "${CHROM_SIZES}" "${SITE_POSITIONS}"; do
  [[ -f "$f" ]] || { echo "Missing reference file: $f"; exit 1; }
done

mkdir -p data/raw data/trimmed results/fastqc_raw results/cutadapt results/hic logs

for SAMPLE in "${!SAMPLES[@]}"; do
  PREFIX="${SAMPLES[$SAMPLE]}"
  log "Processing ${SAMPLE}"

  R1_RAW="data/raw/${SAMPLE}_R1.fastq.gz"
  R2_RAW="data/raw/${SAMPLE}_R2.fastq.gz"

  [[ -f "${R1_RAW}" ]] || wget --no-check-certificate -O "${R1_RAW}" \
    "${BASE_URL}/${PREFIX}_R1_001.fastq.gz"
  [[ -f "${R2_RAW}" ]] || wget --no-check-certificate -O "${R2_RAW}" \
    "${BASE_URL}/${PREFIX}_R2_001.fastq.gz"

  FASTQC_HTML="results/fastqc_raw/${SAMPLE}_R1_fastqc.html"
  if [[ ! -f "${FASTQC_HTML}" ]]; then
    fastqc "${R1_RAW}" "${R2_RAW}" -o results/fastqc_raw --threads "${THREADS}"
  fi

  R1_TRIM="data/trimmed/${SAMPLE}_R1.trimmed.fastq.gz"
  R2_TRIM="data/trimmed/${SAMPLE}_R2.trimmed.fastq.gz"
  CUTADAPT_LOG="results/cutadapt/${SAMPLE}.cutadapt.log"

  if [[ ! -f "${R1_TRIM}" ]]; then
    cutadapt \
      -q "${QUAL_CUTOFF}" \
      -m "${MIN_LEN}" \
      -a "${ADAPTER}" \
      -o "${R1_TRIM}" \
      -p "${R2_TRIM}" \
      "${R1_RAW}" "${R2_RAW}" \
      > "${CUTADAPT_LOG}" 2>&1
  fi

  JUICER_SAMPLE_DIR="data/juicer/${SAMPLE}"
  mkdir -p "${JUICER_SAMPLE_DIR}/fastq"
  ln -sf "$(pwd)/${R1_TRIM}" "${JUICER_SAMPLE_DIR}/fastq/${SAMPLE}_R1.fastq.gz"
  ln -sf "$(pwd)/${R2_TRIM}" "${JUICER_SAMPLE_DIR}/fastq/${SAMPLE}_R2.fastq.gz"

  HIC_OUT="${JUICER_SAMPLE_DIR}/aligned/inter_30.hic"
  if [[ ! -f "${HIC_OUT}" ]]; then
    bash "${JUICER_DIR}/scripts/juicer.sh" \
      -D "${JUICER_DIR}" \
      -d "$(pwd)/${JUICER_SAMPLE_DIR}" \
      -g "${GENOME_SHORT}" \
      -z "$(pwd)/${GENOME_FA}" \
      -p "$(pwd)/${CHROM_SIZES}" \
      -y "$(pwd)/${SITE_POSITIONS}" \
      -s "${ENZYME}" \
      -t "${THREADS}" \
      2>&1 | tee -a "${LOG_DIR}/${SAMPLE}_juicer.log"
  fi

  if [[ -f "${HIC_OUT}" ]]; then
    cp "${HIC_OUT}" "results/hic/${SAMPLE}.inter_30.hic"
  else
    log "WARNING: ${SAMPLE} .hic not found"
  fi
done

log "Done"
ls -lh results/hic/*.hic
