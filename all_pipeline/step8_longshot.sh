#!/usr/bin/env bash
set -euo pipefail

# =======================
# User configs
# =======================
SAMPLE="${1:?Usage: $0 <sample_id> <bam> <ref.fa> <outdir> [threads]}"
BAM="${2:?Usage: $0 <sample_id> <bam> <ref.fa> <outdir> [threads]}"
REF="${3:?Usage: $0 <sample_id> <bam> <ref.fa> <outdir> [threads]}"
OUTDIR="${4:?Usage: $0 <sample_id> <bam> <ref.fa> <outdir> [threads]}"
THREADS="${5:-16}"

# Longshot tuning for ONT + PSMC-friendly conservativeness
MIN_MAPQ=20            # 20-30; use 30 to be conservative for ONT repeats
MIN_ALT_COUNT=4        # 4-6 typical for ONT
MIN_ALT_FRAC=0.2       # 0.2-0.3 to suppress false hets
STRAND_P=0.05          # recommended default for ONT (do not relax)  [3](https://github.com/pjedge/longshot)
QUAL_MIN=30            # filter later
GQ_MIN=30              # filter later

# paths for psmc utils
PSMC_DIR="${PSMC_DIR:-/01_software/Beta-PSMC}"   # change if needed
VCF2SNP_PL="$PSMC_DIR/utils/vcf2snp.pl"
FQ2PSMCFA="$PSMC_DIR/utils/fq2psmcfa"

mkdir -p "$OUTDIR"/{per_chr,tmp,mask,psmc}

# =======================
# Check indexes
# =======================
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"
[[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" ]] || samtools index -@ "$THREADS" "$BAM"

# =======================
# Estimate average depth -> minDP ~ avg/3, maxDP ~ avg*2 (PSMC README suggestion) [1](https://github.com/lh3/psmc)
# =======================
echo "[INFO] Estimating mean depth (this may take a while)..."
MEAN_DP=$(samtools depth -a -q "$MIN_MAPQ" "$BAM" | awk '{sum+=$3; n++} END{if(n>0) printf "%.2f\n", sum/n; else print "0"}')
MIN_DP=$(awk -v d="$MEAN_DP" 'BEGIN{v=int(d/3); if(v<6) v=6; print v}')
MAX_DP=$(awk -v d="$MEAN_DP" 'BEGIN{v=int(d*2); if(v<50) v=50; print v}')

echo "[INFO] MeanDP=$MEAN_DP  MIN_DP=$MIN_DP  MAX_DP=$MAX_DP"

# =======================
# Chromosome list
# =======================
cut -f1 "${REF}.fai" > "$OUTDIR/tmp/contigs.list"

# =======================
# Run Longshot per contig (recommended by author) [3](https://github.com/pjedge/longshot)
# Use -A to auto-set max coverage cutoff (slower but robust) [3](https://github.com/pjedge/longshot)
# =======================
run_longshot_chr () {
  local chr="$1"
  local outvcf="$OUTDIR/per_chr/${SAMPLE}.${chr}.vcf"

longshot \
    --bam "$BAM" \
    --ref "$REF" \
    --out "$outvcf" \
    --sample_id "$SAMPLE" \
    --region "$chr" \
    -A \
    -q "$MIN_MAPQ" \
    -c "$MIN_DP" \
    -e "$MIN_ALT_COUNT" \
    -E "$MIN_ALT_FRAC" \
    -P "$STRAND_P"

  bgzip -f "$outvcf"
  tabix -f -p vcf "${outvcf}.gz"
}
export -f run_longshot_chr
export SAMPLE BAM REF OUTDIR MIN_MAPQ MIN_DP MIN_ALT_COUNT MIN_ALT_FRAC STRAND_P

echo "[INFO] Calling variants with Longshot per contig..."
cat "$OUTDIR/tmp/contigs.list" | xargs -n 1 -P "$THREADS" -I {} bash -c 'run_longshot_chr "$@"' _ {}

# =======================
# Merge per-chr VCFs
# =======================
echo "[INFO] Concatenating per-chr VCFs..."
ls "$OUTDIR/per_chr/${SAMPLE}."*.vcf.gz > "$OUTDIR/tmp/vcfs.list"
bcftools concat -f "$OUTDIR/tmp/vcfs.list" -Oz -o "$OUTDIR/${SAMPLE}.longshot.raw.vcf.gz"
tabix -f -p vcf "$OUTDIR/${SAMPLE}.longshot.raw.vcf.gz"

# =======================
# Normalize + strict filtering for PSMC (keep reliable SNVs)
# - keep biallelic SNP only; PASS; remove dn clusters; QUAL/GQ/DP thresholds
# Longshot marks dense clusters as dn in FILTER field [3](https://github.com/pjedge/longshot)
# =======================
echo "[INFO] Filtering for PSMC..."
bcftools norm -m -any -f "$REF" "$OUTDIR/${SAMPLE}.longshot.raw.vcf.gz" -Ou | \
  bcftools view -m2 -M2 -v snps -Ou | \
  bcftools filter -e 'FILTER~"dn"' -Ou | \
  bcftools view -i "QUAL>=$QUAL_MIN && FORMAT/GQ>=$GQ_MIN && FORMAT/DP>=$MIN_DP && FORMAT/DP<=$MAX_DP" -Oz \
  -o "$OUTDIR/${SAMPLE}.longshot.psmc.filtered.vcf.gz"

tabix -f -p vcf "$OUTDIR/${SAMPLE}.longshot.psmc.filtered.vcf.gz"

# =======================
# Build callable BED mask from depth (mask low/no coverage + ultra-high coverage)
# Important: low coverage should be masked, otherwise consensus may revert-to-ref artifacts [4](https://virological.org/t/issue-with-pipelines-using-bcftools-to-calling-consensus-in-low-coverage-regions/905)
# =======================
echo "[INFO] Building callable regions BED..."
samtools depth -a -q "$MIN_MAPQ" "$BAM" | \
  awk -v min="$MIN_DP" -v max="$MAX_DP" 'BEGIN{OFS="\t"} $3>=min && $3<=max {print $1,$2-1,$2}' \
  > "$OUTDIR/mask/${SAMPLE}.callable.bed"

# Optionally merge adjacent positions to intervals (bedtools recommended if available)
if command -v bedtools >/dev/null 2>&1; then
  sort -k1,1 -k2,2n "$OUTDIR/mask/${SAMPLE}.callable.bed" | bedtools merge -i - \
    > "$OUTDIR/mask/${SAMPLE}.callable.merged.bed"
  CALLABLE_BED="$OUTDIR/mask/${SAMPLE}.callable.merged.bed"
else
  CALLABLE_BED="$OUTDIR/mask/${SAMPLE}.callable.bed"
fi

# =======================
# Generate PSMC input (.psmcfa) using PSMC README dipcall-style route:
# seqtk mutfa ref.fa <(vcf|vcf2snp.pl) | seqtk seq -cM callable.bed | fq2psmcfa -
# [1](https://github.com/lh3/psmc)
# =======================
echo "[INFO] Generating PSMCFA..."
zcat "$OUTDIR/${SAMPLE}.longshot.psmc.filtered.vcf.gz" | \
  "$VCF2SNP_PL" - | \
  seqtk mutfa "$REF" - | \
  seqtk seq -cM "$CALLABLE_BED" -l 80 | \
  "$FQ2PSMCFA" - > "$OUTDIR/psmc/${SAMPLE}.psmcfa"

echo "[DONE] Outputs:"
echo "  - Raw VCF:        $OUTDIR/${SAMPLE}.longshot.raw.vcf.gz"
echo "  - Filtered VCF:   $OUTDIR/${SAMPLE}.longshot.psmc.filtered.vcf.gz"
echo "  - Callable BED:   $CALLABLE_BED"
echo "  - PSMC input:     $OUTDIR/psmc/${SAMPLE}.psmcfa"
