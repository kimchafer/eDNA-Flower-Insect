# 아래 코드로 구성된 파이썬 파일 생성 시 다음의 명령어를 통해 실행 권한을 부여한다: chmod +x filter_m8_by_identity_multi.py [여기서 파일명은 'filter_m8_by_identity_multi.py']
# 스크립트 실행 시 아래와 같이 커맨드를 입력한다
# 1. E2_result.m8 샘플 하나만 수행: python filter_m8_by_identity_multi.py E2_result.m8
# 2. 복수의 샘플에 대해 수행: python filter_m8_by_identity_multi.py E2_result.m8 E3_result.m8 E6_result.m8

import os
import sys

def filter_m8(input_file, output_file):
    """
    Filters .m8 file to retain the highest Identity match for each Query ID.
    If Identity is identical, prioritizes lowest E-value and highest Bit Score.
    """
    query_results = {}

    with open(input_file, 'r') as infile:
        for line in infile:
            cols = line.strip().split('\t')
            query_id, subject_id = cols[0], cols[1]
            identity, e_value, bit_score = float(cols[2]), float(cols[10]), float(cols[11])

            if query_id not in query_results:
                query_results[query_id] = (line, identity, e_value, bit_score)
            else:
                _, max_identity, min_e_value, max_bit_score = query_results[query_id]
                if (identity > max_identity or
                    (identity == max_identity and
                     (e_value < min_e_value or (e_value == min_e_value and bit_score > max_bit_score)))):
                    query_results[query_id] = (line, identity, e_value, bit_score)

    with open(output_file, 'w') as outfile:
        for result in query_results.values():
            outfile.write(result[0])

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python filter_m8_by_identity_multi.py <input_file1.m8> [<input_file2.m8> ...]")
        sys.exit(1)

    for input_file in sys.argv[1:]:
        if not os.path.isfile(input_file):
            print(f"Error: File {input_file} does not exist.")
            continue

        # Create output file name by appending "_filtered" before the extension
        output_file = input_file.replace('.m8', '_filtered.m8')
        print(f"Processing {input_file} -> {output_file}")
        filter_m8(input_file, output_file)
    print("Processing completed.")
