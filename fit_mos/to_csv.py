import re
import sys
import os

def convert_to_csv(input_file):
    # 出力ファイル名は 元の名前_converted.csv
    base, _ = os.path.splitext(input_file)
    output_file = f"{base}_converted.csv"

    try:
        with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        with open(output_file, 'w', encoding='utf-8') as f:
            for line in lines:
                # 1. 前後の余計な空白を削除 (strip)
                # 2. 連続する空白やタブ (\s+) を 1つのカンマに置換
                new_line = re.sub(r'\s+', ',', line.strip())
                f.write(new_line + '\n')
        
        print(f"Success: {output_file} を作成しました。")
        return output_file

    except Exception as e:
        print(f"Error: {e}")
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file1> <input_file2> ...")
    else:
        for file_path in sys.argv[1:]:
            convert_to_csv(file_path)
