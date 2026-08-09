# -*- coding: utf-8 -*-
# 上部にこのマジックコメントを必ず記述してください

# pmosソース接地増幅回路用コード

# データファイル名
data_file = "normal.csv"
data_file_c = "cascode.csv"

# 最小誤差を求めるか
# true: 最小二乗法で最適化を行う
# false: 最適化を行わず、初期値のまま計算
find_min_error = false

# gnuplotの描画をするか。true にすると生成後に gnuplot を呼び出します（PATH に gnuplot が必要）。
run_gnuplot = true

# 係数などの初期値
@k = 0.00034
@vth = -0.7
@R = 2001
@vss = 8.0

# シミュレーション値
k_min = 0.0001
k_max = 0.001
vth_min = -1.2
vth_max = -0.4
step_k = 0.000001
step_vth = 0.001

# 飽和領域のモデル式
def vdp_sat(vgp)
  vdp = @R*0.5*@k*((vgp-@vss-@vth)**2)
  return vdp
end
# 線形領域のモデル式
# 2次方程式 0.5*k*R*vdp^2 - [1 + k*R*(vgp-vth)]*vdp + vss = 0 を解く
def vdp_lin(vgp, vdp_raw)
  a = 0.5 * @k * @R
  b = 1-(@k * @R * (vgp - @vth))
  c = @R*@k*(vgp-0.5*@vss-@vth)*@vss
  
  # 判別式
  discriminant = b**2 - 4*a*c
  
  if discriminant < 0
    # 判別式が負の場合、飽和領域の式を使う
    # （線形領域の境界を超えている可能性がある）
    #puts "Warning: discriminant < 0 at vgp=#{vgp}, using saturation model"
    return vdp_sat(vgp)
  end
  
  # 2つの解のうち、物理的に妥当な方（小さい方）を選択
  vdp1 = (-b - Math.sqrt(discriminant)) / (2*a)
  vdp2 = (-b + Math.sqrt(discriminant)) / (2*a)
  
  # 大きく見ればvdpは単調減少するから、小さい方を選ぶ
  vdp = [vdp1, vdp2].min
  
  return vdp
end
# 線形領域と飽和領域の境界
def vgp_hat(k, vth)
  vgp_hat = ((-1+(1+4*k*@R*@vss)**0.5)/(2*k*@R))+vth
  #puts "vgp_hat=#{vgp_hat} (vgp > vgp_hat: 飽和, vgp <= vgp_hat: 線形)"
  return vgp_hat
end

if find_min_error
  # @kや@vthが動いた時に、最小二乗法で最も近似できているパラメータを決定する。
  for k_test in (k_min..k_max).step(step_k)
    for vth_test in (vth_min..vth_max).step(step_vth)
      @k = k_test
      @vth = vth_test
      error_sum = 0.0
      data = File.new(data_file, "r:UTF-8")
      File.foreach(data) do |line|
        vgp_raw, vdp_raw = line.split.map(&:to_f)
        if vgp_raw <= vgp_hat(@k, @vth) then
          vdp_ideal = vdp_lin(vgp_raw, vdp_raw)
        elsif vgp_raw > @vth
          vdp_ideal = @vss
        else
          vdp_ideal = vdp_sat(vgp_raw)
        end
        error = vdp_ideal - vdp_raw
        error_sum += error**2
      end
      data.close
      #puts "k=#{@k}, vth=#{@vth}, error_sum=#{error_sum}"
      # 誤差最小値の更新
      least_error_sum ||= error_sum
      best_k ||= @k
      best_vth ||= @vth
      if error_sum < least_error_sum
        least_error_sum = error_sum
        best_k = @k
        best_vth = @vth
      end
    end
  end
  puts "best_k=#{best_k}, best_vth=#{best_vth}, least_error_sum=#{least_error_sum}"
  @k = best_k
  @vth = best_vth
end

# 理想曲線を生成
file = File.new("calculated_curve.csv", "w:UTF-8")
data = File.new(data_file, "r:UTF-8")
data_c = File.new(data_file_c, "r:UTF-8")

File.foreach(data) do |line|
  vss_raw, vdp_raw = line.split.map(&:to_f)
  vdp_calculated = vdp_raw / @R
  file.puts "#{vss_raw} #{vdp_calculated}".map(&:to_f)
end
data.close

File.foreach(data) do |line|
  vss_raw_c, vdp_raw_c = line_c.split.map(&:to_f)
end
data_c.close

file.close

output_image = 'vss_id.png'
gp_file = "plot.gp"

gnuplot_script = <<~GP
    set encoding utf8
    set terminal pngcairo font 'Noto Sans,15' size 800,800
    # x軸とy軸のスケールを同じにする
    set size ratio -1
    set datafile separator ' '
    set output '#{output_image}'
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
    # set label "（vss=#{@vss} V, vth=#{format('%.2f', @vth)} V, k=#{(@k*1000000).round} µS/V）" at graph 0.2, -0.2
    # plot '#{data_file}' with linespoints dt 1 pt 7 lc 'black' title '実測値'
    plot '#{data_file}' with linespoints dt 1 pt 7 lc 'black' title 'PMOS1つ', '#{data_file_c}' with linespoints dt 2 lc 'black' pt 9 title 'PMOSカスコード接続'
GP

# UTF-8 で書き出す
File.open(gp_file, "w:UTF-8") do |f|
    f.write(gnuplot_script)
end

puts "wrote #{gp_file} -> #{output_image}"

if run_gnuplot
    # 実行したい場合は PATH に gnuplot が通っていることを確認してください
    system("gnuplot", gp_file)
end