# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file
set grid

# --- グラフ装飾 ---
set xlabel "X-axis"
set ylabel "Y-axis"
set title sprintf("Data Plot: %s", file_A)

# --- プロットの実行 ---
# using 1:2 で1列目を横軸、2列目を縦軸に指定
# with linespoints で点と線の両方を表示（お好みで with points や with lines に変更可）
plot file_A using 1:2 title "Measured Data" with linespoints pt 7 lc rgb "dark-violet"

print "--- Simple Plot Generated ---"
print "File: " . file_A
print "Output: " . out_file
