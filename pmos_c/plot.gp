set encoding utf8
set terminal pngcairo font 'Noto Sans,15' size 800,800
# x軸とy軸のスケールを同じにする
set size ratio -1
set datafile separator ' '
set output 'normal.png'
unset title
set grid
set xlabel "電源電圧vss [V]"
set ylabel "出力電圧vdp [V]"
set xrange [4:8]
set yrange [0:3]
set xtics 1
set mxtics 5
set ytics 1
set mytics 5
set bmargin 6
# set label "（vss=8.0 V, vth=-0.70 V, k=340 µS/V）" at graph 0.2, -0.2
# plot 'normal.csv' with linespoints dt 1 pt 7 lc 'black' title '実測値'
plot 'normal.csv' with linespoints dt 1 pt 7 lc 'black' title 'PMOS1つ', 'cascode.csv' with linespoints dt 2 lc 'black' pt 9 title 'PMOSカスコード接続'
