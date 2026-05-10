# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file
set grid

# --- グラフ装飾 ---
set xlabel "X-axis"
set ylabel "Y/Z-axis"
set title "Comparison of Two Data Sets"

# --- プロットの実行 ---
# file_B が空文字列でないか確認して分岐
if (exists("file_B") && strlen(file_B) > 0) {
    plot file_A using 1:2 title sprintf("File A: %s", file_A) with points pt 7 lc rgb "dark-spring-green", \
         file_B using 1:2 title sprintf("File B: %s", file_B) with lines lw 2 lc rgb "web-blue"
} else {
    print "Warning: file_B is not specified. Plotting file_A only."
    plot file_A using 1:2 title sprintf("File A: %s", file_A) with points pt 7 lc rgb "dark-spring-green"
}
