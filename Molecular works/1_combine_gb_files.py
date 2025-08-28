import os

# 입력 폴더와 출력 파일 경로 설정
input_folder = "/storage2/jihoonkim/eDNA2/DB/gb_files"  # 여러 .gb 파일이 저장된 폴더 경로
output_file = "COI_all_Insecta.gb"         # 합쳐진 파일의 출력 경로

# 파일 합치기
with open(output_file, 'w') as outfile:
    for filename in os.listdir(input_folder):
        if filename.endswith(".gb"):  # .gb 확장자 확인
            file_path = os.path.join(input_folder, filename)
            with open(file_path, 'r') as infile:
                content = infile.read()
                # 각 파일 내용을 추가
                outfile.write(content)
                # 구분자 추가 (중복 방지용으로 "//" 추가 확인)
                if not content.strip().endswith("//"):
                    outfile.write("\n//\n")

print(f"모든 .gb 파일이 {output_file}로 합쳐졌습니다.")
