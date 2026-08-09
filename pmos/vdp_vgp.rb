# -*- coding: utf-8 -*-
# 上部にこのマジックコメントを必ず記述してください

# pmosソース接地増幅回路用コード

# データファイル名
data_file = "vgp_vdp_0106.csv"

# 最小誤差を求めるか
# true: 最小二乗法で最適化を行う
# false: 最適化を行わず、初期値のまま計算
find_min_error = false

# gnuplotの描画をするか。true にすると生成後に gnuplot を呼び出します（PATH に gnuplot が必要）。
run_gnuplot = true

# 係数などの初期値
@k = 0.000195
@vth = -0.4
@R = 2001
@vss = 8.0

# シミュレーション値
k_min = 0.00001
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
def vdp_lin(vgp)
  a = 0.5 * @k * @R
  b = 1-(@k * @R * (vgp - @vth))
  c = @R*@k*(vgp-0.5*@vss-@vth)*@vss
  
  # 判別式
  discriminant = b**2 - 4*a*c
  
  if discriminant < 0
    # 判別式が負の場合、飽和領域の式を使う
    # （線形領域の境界を超えている可能性がある）
    puts "Warning: discriminant < 0 at vgp=#{vgp}, using saturation model"
    return vdp_sat(vgp)
  end
  
  # 2つの解のうち、大きい方でないとvdpが負になる
  vdp = (-b + (discriminant)**0.5) / (2*a)
  
  return vdp
end

# 線形領域と飽和領域の境界
def vgp_hat(k, vth)
  # 図を描けば、飽和領域が下に凸、条件式の境界の方程式が単調増加であるから、vgp_hatは1つに決まる。
  vgp_hat = ((1+k*@R*@vss-(1+2*k*@R*@vss)**0.5)/(k*@R))+vth
  # puts "vgp_hat=#{vgp_hat} (vgp <= vgp_hat: 飽和, vgp < vgp_hat: 線形)"
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
         vdp_ideal = vdp_lin(vgp_raw)
        elsif vgp_raw > @vss + @vth
         vdp_ideal = 0.0
        else
         vdp_ideal = vdp_sat(vgp_raw)
        end

        # if vgp_raw <= 3.0 then
        #  vdp_ideal = vdp_lin(vgp_raw)
        # elsif vgp_raw > 7.0
        #  vdp_ideal = 0.0
        # else
        #  vdp_ideal = vdp_sat(vgp_raw)
        # end
        
        # if vgp_raw <= 3.0 then
        #   vdp_ideal = vdp_sat(vgp_raw)
        # elsif vgp_raw > 7.0
        #   vdp_ideal = 0.0
        # else
        #   vdp_ideal = vdp_lin(vgp_raw)
        # end
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

File.foreach(data) do |line|
  vgp_raw, vdp_raw = line.split.map(&:to_f)
  if vgp_raw <= vgp_hat(@k, @vth) then
    vdp_ideal = vdp_lin(vgp_raw)
    log = "lin"
  elsif vgp_raw > @vss + @vth
    vdp_ideal = 0.0
    log = "cutoff"
  else
    vdp_ideal = vdp_sat(vgp_raw)
    log = "sat"
  end

  # if vgp_raw <= 3.0 then
  #   vdp_ideal = vdp_lin(vgp_raw)
  #   log = "lin"
  # elsif vgp_raw > 7.0
  #   vdp_ideal = 0.0
  #   log = "cutoff"
  # else
  #   vdp_ideal = vdp_sat(vgp_raw)
  #   log = "sat"
  # end
  
  # if vgp_raw <= 3.0 then
  #   vdp_ideal = vdp_sat(vgp_raw)
  #   log = "sat"
  # elsif vgp_raw > 7.0
  #   vdp_ideal = 0.0
  #   log = "cutoff"
  # else
  #   vdp_ideal = vdp_lin(vgp_raw)
  #   log = "lin"
  # end
  file.write(vgp_raw)
  file.write(" ")
  file.write(vdp_ideal)
  file.write(" ")
  file.puts(log)
end
data.close
file.close

puts "vgp_hat=#{vgp_hat(@k, @vth)}"

output_image = data_file.sub(/\.csv$/, '.png')
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
    set xlabel "入力電圧vgp [V]"
    set ylabel "出力電圧vdp [V]"
    #set xrange [0:8]
    #set yrange [0:8]
    set xtics 1
    set mxtics 5
    set ytics 1
    set mytics 5
    set bmargin 6
    set label "（vss=#{@vss} V, vth=#{format('%.2f', @vth)} V, k=#{(@k*1000000).round} µS/V, vgp hat=#{format('%.2f', vgp_hat(@k, @vth))} V）" at graph 0.1, -0.2
    #plot '#{data_file}' with linespoints dt 2 lc 'black' title 'plot'
    plot '#{data_file}' lc 'black' pt 7 title '実測値', 'calculated_curve.csv' with linespoints dt 2 lc 'black' pt 9 title '計算値'
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