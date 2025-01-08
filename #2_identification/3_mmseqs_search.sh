#!/bin/bash
#SBATCH --job-name=mmseqs_search
#SBATCH --output=mmseqs_search_%j.out
#SBATCH --error=mmseqs_search_%j.err
#SBATCH --time=UNLIMITED
#SBATCH -p Node7
#SBATCH -n 16

# 본 스크립트는 각 샘플마다 개별(병렬)적으로 분석을 수행하기 위해 단일 샘플을 위해 작성되었다.
# 스크립트 실행 시 아래와 같이 커맨드를 입력한다
# sbatch 3_mmseqs_search.sh E2 (E2는 raw data의 샘플명에 해당한다. raw data의 full name은 시퀀싱 기관에 따라 상이할 수 있기에 아래 1번 단계의 양식을 유의하라)

# 입력 변수
sample=$1
if [ -z "$sample" ]; then
  echo "Error: No sample name provided. Run as: sbatch script_name.sh <sample>"
  exit 1
fi

# 작업 시작 시간 기록
start_time=$(date +%s)

# 디렉토리 설정
work_dir=/storage2/jihoonkim/eDNA/analysis/2_20250106  # 사전에 만들어진 워킹 디렉토리
raw_data_dir=${work_dir}/raw_data                      # 사전에 만들어진 raw data 디렉토리
processed_data_dir=${work_dir}/processed_data
results_dir=${work_dir}/results
tmp_dir=${work_dir}/tmp
ref_db=/storage2/jihoonkim/eDNA/mmseqs/mmseqs_cox1/cox1_refDB  # Reference DB (MMseqs2 형식)

mkdir -p $processed_data_dir $results_dir $tmp_dir

# 1. Paired-End 데이터 병합
echo "==> Merging paired-end reads for sample: $sample"  # 아래 .fastq.gz 형식의 raw data는 시퀀싱 기관으로부터 받은 것이며, 기관마다 해당 이름들의 양식이 상이할 수 있기에, 경우에 따라 아래 코드를 수정할 필요가 있다.현재 기준 SNU NICEM
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
