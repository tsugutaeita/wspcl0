set encoding utf8
set terminal pngcairo font 'Noto Sans,15' size 800,800
# x軸とy軸のスケールを同じにする
set size ratio -1
set datafile separator ' '
set output 'vgs_vds_1221.png'
unset title
set grid
set xlabel "ゲート・ソース間電圧VGS [V]"
set ylabel "ドレイン・ソース間電圧VDS [V]"
set xrange [0:8]
set yrange [0:7]
set xtics 1
set mxtics 5
set ytics 1
set mytics 5
set bmargin 6
set label "VGS-VDS特性曲線（vdd=6.0 V, vth=0.516 V, k=338 µS/V）" at graph 0.05, -0.2
#plot 'vgs_vds_1221.csv' with linespoints dt 2 lc 'black' title 'plot'
plot 'vgs_vds_1221.csv' lc 'black' pt 7 title '実測値', 'calculated_curve.csv' with linespoints dt 2 lc 'black' pt 9 title '近似値'
