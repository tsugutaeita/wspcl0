# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file
set grid lc rgb "black" lw 3
set border lc rgb "black" lw 3

# --- グラフ装飾 ---
set xlabel "X-axis"
set ylabel "Y-axis"
set title sprintf("Data Plot: %s", file_a) tc rgb "black"

# テキストボックスのスタイル設定（背景を白、枠線なし）
set style textbox opaque fc rgb "white" noborder

# 凡例（key）の位置を調整（ラベルと重なる場合は右下などに移動）
set key bottom right opaque

# セパレータ設定（前回の解決策を反映）
if (exists("sep") && sep eq ",") {
    set datafile separator ","
} else {
    set datafile separator whitespace
}

# --- プロットの実行 ---
# using 1:2 で1列目を横軸、2列目を縦軸に指定
# with linespoints で点と線の両方を表示（お好みで with points や with lines に変更可）
plot file_a using 1:2 title "Measured Data" with points pt 7 lc rgb "black"

print "--- Simple Plot Generated ---"
print "File: " . file_a
print "Output: " . out_file
