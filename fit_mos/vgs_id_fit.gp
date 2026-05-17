# --- 設定 ---
reset
set terminal pngcairo size w_size, h_size
set output out_file
set grid lc rgb "black" lw 3
set border lc rgb "black" lw 3

# --- フィッティング定義 ---
# モデル式を文字列としても定義しておく（図中出力用）
model_formula = "I_D = (k/2) * (V_{GS} - V_{th})^2"
f(x) = (x > Vth) ? (k/2.0) * (x - Vth)**2 : 0

# Pythonから渡された初期値をセット
k = k_i
Vth = vth_i

# セパレータ設定（前回の解決策を反映）
if (exists("sep") && sep eq ",") {
    set datafile separator ","
} else {
    set datafile separator whitespace
}

# テキストボックスのスタイル設定（背景を白、枠線なし）
set style textbox opaque fc rgb "white" noborder

# label の末尾に "boxed" を追加
set label 1 sprintf("Model: %s\nV_{th} = %.3f [V]\nk = %.1f [uS/V]", model_formula, Vth, k_uSV) \
          at graph 0.05, 0.9 \
          font ",10" tc rgb "red" boxed

# フィッティング実行
fit f(x) file_a using 1:2 via k, Vth

# --- 図中ラベルの設定 ---
# graph 座標系（0.0〜1.0）を使用して左上（0.05, 0.9）あたりに表示
# \n は改行コードです
set label 1 sprintf("Model: %s\nV_{th} = %.3f [V]\nk = %.2e [A/V^2]", model_formula, Vth, k) \
          at graph 0.05, 0.9 \
          font ",10" tc rgb "black"

# --- プロット ---
set xlabel "V_{GS} [V]"
set ylabel "I_D [A]"
set title sprintf("MOSFET Parameter Fitting\nSource: %s", file_a) noenhanced

# 凡例（key）の位置を調整（ラベルと重なる場合は右下などに移動）
set key bottom right opaque

plot file_A using 1:2 title "Measured Data" with points pt 7 lc rgb "black", \
     f(x) title "Fitted Curve" with lines lw 2 lc rgb "black" dt 1

# コンソールへの結果表示
print "========================================="
print "  Fitting Results"
print "========================================="
print sprintf("  Vth : %f [V]", Vth)
print sprintf("  k   : %e [A/V^2]", k)
print "========================================="
