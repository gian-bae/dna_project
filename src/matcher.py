#src/matcher.py

def load_fasta(file_path):
    """
    FASTA 파일을 읽어서 헤더(Header)와 시퀀스(Sequence)를 딕셔너리로 파싱하는 함수입니다.
    """
    sequences = {}
    current_header = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue #빈 줄은 무시

                if line.startswith('>'):
                    current_header = line[1:] # '>' 기호 제거하고 키값으로 설정
                    sequences[current_header] = ""
                else:
                    if current_header:
                        sequences[current_header] += line.upper() #대문자로 통일
                        
        return sequences
                    
    except FileNotFoundError:
        print(f"[Error] 파일을 찾을 수 없습니다: {file_path}")
        return {}

#테스트 실행 코드 (이 파일을 직접 실행했을 때만 동작)
if __name__ == "__main__":
    #상위 폴더의 data 디렉토리 경로 지정
    fasta_data = load_fasta ("data/sample.fasta")

    print("== FASTA Parsing Result ===")
    for header, seq in fasta_data.items():
        print(f"Header: {header}")
        print(f"Length: {len(seq)} bp")
        print(f"Sequence: {seq}\n")
