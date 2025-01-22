#!/bin/bash
#SBATCH --job-name=full_analysis_pipeline    # Job name
#SBATCH --output=full_analysis_%j.out       # Standard output
#SBATCH --error=full_analysis_%j.err        # Error output
#SBATCH --time=UNLIMITED                    # No time limit
#SBATCH -p Node7                            # Partition
#SBATCH -n 32                               # Number of CPUs

# 샘플 이름 설정
sample=$1
if [ -z "$sample" ]; then
  echo "Error: No sample name provided. Run as: sbatch script_name.sh <sample>"
  exit 1
fi

# 작업 시작 시간 기록
start_time=$(date +%s)

# 디렉토리 설정
work_dir=/storage2/jihoonkim/eDNA2/analysis/2_20250120
raw_data_dir=${work_dir}/raw_data
processed_data_dir=${work_dir}/processed_data
results_dir=${work_dir}/results
tmp_dir=${work_dir}/tmp
ref_db=/storage2/jihoonkim/eDNA2/DB/cox1.refDB  # Reference DB

mkdir -p $processed_data_dir $results_dir $tmp_dir

# 1. Paired-End 데이터 병합
echo "==> Merging paired-end reads for sample: $sample"
R1=$(find $raw_data_dir -name "${sample}_S*_R1_001.fastq.gz")
R2=$(find $raw_data_dir -name "${sample}_S*_R2_001.fastq.gz")
merged_prefix=${processed_data_dir}/${sample}_merged

pear -f $R1 -r $R2 -o $merged_prefix -j 8

merged_file=${merged_prefix}.assembled.fastq
if [ ! -f "$merged_file" ]; then
  echo "Error: Merged file not found for sample $sample."
  exit 1
fi

# 2. FASTQ -> FASTA 변환
echo "==> Converting FASTQ to FASTA for sample: $sample"
fasta_file=${processed_data_dir}/${sample}_merged.fasta
seqtk seq -A $merged_file > $fasta_file

# 3. MMseqs2 데이터베이스 생성
echo "==> Creating MMseqs2 database for sample: $sample"
query_db=${processed_data_dir}/${sample}_queryDB
mmseqs createdb $fasta_file $query_db

# 4. Taxonomy 결과 경로 설정
taxonomy_result=${results_dir}/${sample}_taxonomy_result
tsv_result=${results_dir}/${sample}_taxonomy_result.tsv
krona_report=${results_dir}/${sample}_krona_report.html
taxonomy_tmp_dir=${tmp_dir}/${sample}_taxonomy_tmp
mkdir -p $taxonomy_tmp_dir

# 5. MMseqs2 Taxonomy 분석 수행
echo "==> Running taxonomy analysis for sample: $sample"
mmseqs taxonomy $query_db $ref_db $taxonomy_result $taxonomy_tmp_dir \
  --tax-lineage 1 --threads 32 --search-type 3 --orf-filter 0 --min-seq-id 0.95 --alignment-mode 4
if [ $? -ne 0 ]; then
  echo "Error: Taxonomy analysis failed for sample $sample."
  exit 1
fi

# 6. Taxonomy 결과를 TSV로 변환
echo "==> Converting taxonomy results to TSV for sample: $sample"
mmseqs createtsv $query_db $taxonomy_result $tsv_result --threads 16
if [ $? -ne 0 ]; then
  echo "Error: Failed to convert taxonomy results to TSV for sample $sample."
  exit 1
fi

# 7. Krona HTML 보고서 생성
echo "==> Generating Krona HTML report for sample: $sample"
mmseqs taxonomyreport $ref_db $taxonomy_result $krona_report --threads 16 --report-mode 1
if [ $? -ne 0 ]; then
  echo "Error: Failed to generate Krona report for sample $sample."
  exit 1
fi

# 작업 종료 시간 기록 및 소요 시간 계산
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "==> Full analysis pipeline completed for sample: $sample"
echo "Results summary:"
echo "  - Taxonomy DB: $taxonomy_result"
echo "  - TSV Result: $tsv_result"
echo "  - Krona HTML Report: $krona_report"
echo "Total elapsed time: $((elapsed_time / 60)) minutes and $((elapsed_time % 60)) seconds"
