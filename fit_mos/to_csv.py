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
        # 引数で渡されたファイルパスのリストを取得
        file_paths = sys.argv[1:]
        
        # 変換処理に回すファイルを格納するリスト
        targets_to_convert = []

        # 1. 先に末尾が _converted.csv のファイルだけを処理（削除）
        for file_path in file_paths:
            if file_path.endswith('_converted.csv'):
                if os.path.exists(file_path):
                    try:
                        os.remove(file_path)
                        print(f"Removed: {file_path} (古い変換済ファイルを削除しました)")
                    except Exception as e:
                        print(f"Error removing {file_path}: {e}")
            else:
                # _converted.csv 以外のファイルは後で変換するためキープ
                targets_to_convert.append(file_path)

        # 2. 残ったファイルに対してのみ、順次変換処理を実行
        for file_path in targets_to_convert:
            if os.path.exists(file_path):
                convert_to_csv(file_path)
            else:
                print(f"Warning: {file_path} が見つかりません。スキップします。")