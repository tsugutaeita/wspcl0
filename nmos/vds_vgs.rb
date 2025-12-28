# -*- coding: utf-8 -*-
# 上部にこのマジックコメントを必ず記述してください

# nmosソース接地増幅回路用コード

# データファイル名
data_file = "vgs_vds_1221.csv"

# 最小誤差を求めるか
# true: 最小二乗法で最適化を行う
# false: 最適化を行わず、初期値のまま計算
find_min_error = true

# 係数などの初期値
# 0.5*unCox*W/Lをkとして定義
@k = 0.00034
@vth = 0.5
@R = 2011
@vdd = 6.0

# シミュレーション値
k_min = 0.0001
k_max = 0.001
vth_min = 0.4
vth_max = 1.0
step_k = 0.000001
step_vth = 0.001

# 飽和領域のモデル式
def vds_sat(vgs)
  vds = @vdd-@R*@k*((vgs-@vth)**2)
  return vds
end
# 線形領域のモデル式
# 2次方程式 k*R*vds^2 - (1 + 2*k*R*(vgs-vth))*vds + vdd = 0 を解く
def vds_lin(vgs, vds_raw)
  a = @k*@R
  b = -(1 + 2 * @k * @R * (vgs - @vth))
  c = @vdd
  
  # 判別式
  discriminant = b**2 - 4*a*c
  
  if discriminant < 0
    # 判別式が負の場合、飽和領域の式を使う
    # （線形領域の境界を超えている可能性がある）
    #puts "Warning: discriminant < 0 at vgs=#{vgs}, using saturation model"
    return vds_sat(vgs)
  end
  
  # 解析の結果、2つの解のうち小さい方が妥当と判明した。（Vds < Vgs-Vth）
  vds = (-b - Math.sqrt(discriminant)) / (2*a)
  
  return vds
end
# 線形領域と飽和領域の境界
def vgs_hat(k, vth)
  vgs_hat = ((-1+(1+4*k*@R*@vdd)**0.5)/(2*k*@R))+vth
  #puts "vgs_hat=#{vgs_hat} (VGS < vgs_hat: 飽和, VGS >= vgs_hat: 線形)"
  return vgs_hat
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
        vgs_raw, vds_raw = line.split.map(&:to_f)
        if vgs_raw >= vgs_hat(@k, @vth) then
          vds_ideal = vds_lin(vgs_raw, vds_raw)
        elsif vgs_raw < @vth
          vds_ideal = @vdd
        else
          vds_ideal = vds_sat(vgs_raw)
        end
        error = vds_ideal - vds_raw
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

File.foreach(data) do |line|
  vgs_raw, vds_raw = line.split.map(&:to_f)
  if vgs_raw >= vgs_hat(@k, @vth) then
    vds_ideal = vds_lin(vgs_raw, vds_raw)
  #  log = "lin"
  elsif vgs_raw < @vth
    vds_ideal = @vdd
  else
    vds_ideal = vds_sat(vgs_raw)
  #  log = "sat"
  end
  file.write(vgs_raw)
  file.write(" ")
  file.puts(vds_ideal)
  #file.write(" ")
  #file.puts(log)
end
data.close
file.close

# 実行するかどうか。true にすると生成後に gnuplot を呼び出します（PATH に gnuplot が必要）。
run_gnuplot = true

output_image = 'vgs_vds_1221.png'
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
    set xlabel "ゲート・ソース間電圧VGS [V]"
    set ylabel "ドレイン・ソース間電圧VDS [V]"
    set xrange [0:8]
    set yrange [0:7]
    set xtics 1
    set mxtics 5
    set ytics 1
    set mytics 5
    set bmargin 6
    set label "VGS-VDS特性曲線（vdd=#{@vdd} V, vth=#{@vth} V, k=#{@k} S/V）" at graph 0.05, -0.2
    #plot '#{data_file}' with linespoints dt 2 lc 'black' title 'plot'
    plot '#{data_file}' lc 'black' pt 7 title '実測値', 'calculated_curve.csv' with linespoints dt 2 lc 'black' pt 9 title '近似値'
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