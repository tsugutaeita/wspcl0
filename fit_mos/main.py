import os
import subprocess
import glob

def run_interpreter():
    # 1. 実行する .gp ファイルの選択
    gp_files = glob.glob("*.gp")
    if not gp_files:
        print("Error: .gpファイルが見つかりません。")
        return

    print("--- [1/3] 実行するGnuplotスクリプトを選択 ---")
    for i, f in enumerate(gp_files):
        print(f"[{i}] {f}")
    
    try:
        gp_idx = int(input("番号入力 > "))
        selected_gp = gp_files[gp_idx]
    except: return

    # 2. 解析対象データの選択
    csv_files = glob.glob("csv/*.csv")
    if not csv_files:
        print("Error: csv/ ディレクトリにファイルがありません。")
        return
    
    print("\n--- [2/3] データファイルA を選択 ---")
    for i, f in enumerate(csv_files):
        print(f"[{i}] {f}")
    try:
        idx_a = int(input("番号入力 > "))
        file_a = csv_files[idx_a]
    except: return

    print("\n--- [3/3] データファイルB を選択 (不要な場合はEnter) ---")
    for i, f in enumerate(csv_files):
        print(f"[{i}] {f}")
    
    choice_b = input("番号入力 > ")
    file_b = csv_files[int(choice_b)] if choice_b.strip() else ""

    # 3. パラメータの入力
    print("\n--- オプション設定 (Enterでデフォルト値) ---")
    
    # 区切り文字の設定
    print("データファイルA用の区切り文字を選択: [0] カンマ(,)  [1] タブ(\\t)  [2] スペース( )")
    sep_choice = input("番号入力 [default: 0] > ") or "0"
    sep_map = {"0": ",", "1": "\\t", "2": " "}
    sep_a = sep_map.get(sep_choice, ",")
    
    # --- 修正箇所：file_b がある場合のみ質問する ---
    sep_b = "" 
    if file_b:
        print("データファイルB用の区切り文字を選択: [0] カンマ(,)  [1] タブ(\\t)  [2] スペース( )")
        sep_choice_b = input("番号入力 [default: 0] > ") or "0"
        sep_b = sep_map.get(sep_choice_b, ",")
    # ----------------------------------------------

    base_name = os.path.splitext(os.path.basename(file_a))[0]
    out_file = input(f"出力ファイル名 [default: {base_name}_out.png] > ") or f"{base_name}_out.png"
    vth_init = input("Vth 初期値 [default: 0.5] > ") or "0.5"
    k_init = input("k 初期値 [default: 0.001] > ") or "0.001"
    width = input("画像幅 [default: 800] > ") or "800"
    height = input("画像高さ [default: 600] > ") or "600"

    # 4. Gnuplot引数の組み立て
    # sepを変数として追加。文字列なのでクォートで囲む
    e_string = (
        f"file_a='{file_a}'; "
        f"file_b='{file_b}'; "
        f"out_file='{out_file}'; "
        f"sep_a='{sep_a}'; "
        f"sep_b='{sep_b}'; "
        f"w_size={width}; "
        f"h_size={height}; "
        f"vth_i={vth_init}; "
        f"k_i={k_init}"
    )

    # 5. 実行
    print(f"\n実行コマンド: gnuplot -e \"{e_string}\" {selected_gp}")
    
    try:
        result = subprocess.run(
            ["gnuplot", "-e", e_string, selected_gp],
            capture_output=True, text=True, check=True
        )
        print("\n--- Gnuplot 実行結果 ---")
        if result.stdout: print(result.stdout)
        print(f"完了: {out_file} を出力しました。")
        
    except subprocess.CalledProcessError as e:
        print(f"\n--- Gnuplot エラー発生 ---")
        print(e.stderr)

if __name__ == "__main__":
    if not os.path.exists("csv"): os.makedirs("csv")
    run_interpreter()
