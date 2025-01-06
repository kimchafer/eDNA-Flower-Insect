#!/bin/bash
#SBATCH --job-name=mmseqs_search
#SBATCH --output=mmseqs_search_%j.out
#SBATCH --error=mmseqs_search_%j.err
#SBATCH --mem=64G
#SBATCH --time=UNLIMITED
#SBATCH --cpus-per-task=8
#SBATCH --partition=Node3

# 아래 코드로 구성된 파이썬 파일 생성 시 다음의 명령어를 통해 실행 권한을 부여한다: chmod +x filter_m8_by_identity_multi.py [여기서 파일명은 'filter_m8_by_identity_multi.py']
# 스크립트 실행 시 아래와 같이 커맨드를 입력한다
# 1. E2_result.m8 샘플 하나만 수행: python filter_m8_by_identity_multi.py E2_result.m8
# 2. 복수의 샘플에 대해 수행: python filter_m8_by_identity_multi.py E2_result.m8 E3_result.m8 E6_result.m8

# 입력 변수
sample=$1
if [ -z "$sample" ]; then
  echo "Error: No sample name provided. Run as: sbatch script_name.sh <sample>"
  exit 1
fi

# 작업 시작 시간 기록
start_time=$(date +%s)

# 디렉토리 설정
work_dir=/storage2/jihoonkim/eDNA/analysis/1_20250104
raw_data_dir=${work_dir}/raw_data
processed_data_dir=${work_dir}/processed_data
results_dir=${work_dir}/results
tmp_dir=${work_dir}/tmp
ref_db=/storage2/jihoonkim/eDNA/data/co1/cox1_refDB  # Reference DB (MMseqs2 형식)

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

# 4. MMseqs2 검색 수행
echo "==> Running MMseqs2 search for sample: $sample"
result_db=${results_dir}/${sample}_resultDB
mmseqs search $query_db $ref_db $result_db $tmp_dir --threads 8 -e 1e-5 --search-type 3 --min-seq-id 0.9

# 5. 결과 변환
echo "==> Converting results to .m8 format for sample: $sample"
output_m8=${results_dir}/${sample}_result.m8
mmseqs convertalis $query_db $ref_db $result_db $output_m8

echo "==> Analysis completed for sample: $sample"

# 작업 종료 시간 기록 및 총 소요 시간 계산
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Total elapsed time: $((elapsed_time / 60)) minutes and $((elapsed_time % 60)) seconds"
