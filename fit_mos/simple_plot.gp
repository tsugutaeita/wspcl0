# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file
set grid lc rgb "black" lw 3
set border lc rgb "black" lw 3

# --- グラフ装飾 ---
# 横軸・縦軸ともに2列目・3列目のデータ(Volt)に合わせてラベルを設定
set xlabel "Voltage 1 [V]"
set ylabel "Voltage 2 [V]"
set title sprintf("Data Plot: %s", file_a) tc rgb "black"

# テキストボックスのスタイル設定（背景を白、枠線なし）
set style textbox opaque fc rgb "white" noborder

# 凡例（key）の位置を調整
set key bottom right opaque

# セパレータ設定（CSVデータに対応するためデフォルトまたは指定時にカンマ区切りに設定）
if (exists("sep") && sep eq ",") {
    set datafile separator ","
} else {
    set datafile separator ","
}

# --- プロットの実行 ---
# データヘッダー（最初の2行）をスキップし、using 2:3 で2列目を横軸、3列目を縦軸に指定
plot file_a skip 2 using 2:3 title "Measured Data" with points pt 7 lc rgb "black"

print "--- Simple Plot Generated ---"
print "File: " . file_a
print "Output: " . out_file
