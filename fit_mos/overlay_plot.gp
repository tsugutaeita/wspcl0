# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file
set grid

# 全ファイルをカンマ区切りとして扱う設定
set datafile separator ","

# --- グラフ装飾 ---
set xlabel "X-axis"
set ylabel "Y/Z-axis"
set title sprintf("Data Plot: %s %s", file_a, file_b) tc rgb "black"

# テキストボックスと凡例の設定
set style textbox opaque fc rgb "white" noborder
set key top left opaque

# --- プロットの実行 ---
# file_b が存在し、かつ空文字列でない場合のみ両方をプロット
if (exists("file_b") && strlen(file_b) > 0) {
    plot file_a using 1:2 title sprintf("File A: %s", file_a) with points pt 7 lc rgb "red", \
         file_b using 1:2 title sprintf("File B: %s", file_b) with points pt 7 lc rgb "blue"
} else {
    # file_b が指定されていない場合は file_a のみプロット
    plot file_a using 1:2 title sprintf("File A: %s", file_a) with points pt 7 lc rgb "red"
}
