set encoding utf8
set terminal pngcairo font 'Noto Sans,15' size 800,800
# x軸とy軸のスケールを同じにする
set size ratio -1
set datafile separator ' '
set output 'vgp_vdp_0106.png'
unset title
set grid
set xlabel "入力電圧vgp [V]"
set ylabel "出力電圧vdp [V]"
#set xrange [0:8]
#set yrange [0:8]
set xtics 1
set mxtics 5
set ytics 1
set mytics 5
set bmargin 6
set label "（vss=8.0 V, vth=-0.40 V, k=195 µS/V, vgp hat=3.27 V）" at graph 0.1, -0.2
#plot 'vgp_vdp_0106.csv' with linespoints dt 2 lc 'black' title 'plot'
plot 'vgp_vdp_0106.csv' lc 'black' pt 7 title '実測値', 'calculated_curve.csv' with linespoints dt 2 lc 'black' pt 9 title '計算値'
