# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file

# 参照元の美しい黒枠と太いグリッド線を採用
set grid lc rgb "black" lw 3
set border lc rgb "black" lw 3

# カンマ区切り設定
set datafile separator ","

# --- グラフ装飾 ---
set xlabel "V_{GS} [V]"
set ylabel "I_D [A]"
#set yrange [-0.008:0.008]

# file_b の有無でタイトルを動的に切り替え
if (exists("file_b") && strlen(file_b) > 0) {
    set title sprintf("MOSFET Data Comparison\nSources: %s & %s", file_a, file_b) noenhanced
} else {
    set title sprintf("MOSFET Data Plot\nSource: %s", file_a) noenhanced
}

# テキストボックスと凡例の設定
set style textbox opaque fc rgb "white" noborder
set key top right opaque

# --- プロットの実行 ---
# file_b を先に、file_a を後に書くことで、file_a が一番上に重なります
if (exists("file_b") && strlen(file_b) > 0) {
    plot file_b using 1:2 title "TL082ACPパラメータでのシミュレーション" with points pt 7 lc rgb "blue", \
         file_a using 1:2 title "モデル化したオペアンプ" with points pt 7 lc rgb "red"
} else {
    # file_b が指定されていない場合は file_a のみプロット
    plot file_a using 1:2 title "Measured Data A" with points pt 7 lc rgb "red"
}
