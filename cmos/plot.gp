set encoding utf8
set terminal pngcairo font 'Noto Sans,15' size 800,800
# x軸とy軸のスケールを同じにする
set size ratio -1
set datafile separator ' '
set output 'vin_vout_0112.png'
unset title
set grid
set xlabel "入力電圧vin [V]"
set ylabel "出力電圧vout [V]"
set xrange [0:8]
set yrange [0:8]
set xtics 1
set mxtics 5
set ytics 1
set mytics 5
set bmargin 6
# set label "（vss=8.0 V, vth=-0.70 V, k=340 µS/V）" at graph 0.2, -0.2
plot 'vin_vout_0112.csv' with linespoints dt 1 pt 7 lc 'black' title '実測値'
# plot 'vin_vout_0112.csv' lc 'black' pt 7 title '実測値', 'calculated_curve.csv' with linespoints dt 2 lc 'black' pt 9 title '計算値'
