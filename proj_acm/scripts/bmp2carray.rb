#!/usr/bin/env ruby

def convert_bmp_to_c_array(file_path)
  File.open(file_path, "rb") do |f|
    # BMPヘッダ(14バイト)とDIBヘッダ(40バイト)を読み込む
    header = f.read(54)
    if header[0, 2] != "BM"
      abort "エラー: 有効なBMPファイルではありません。"
    end

    # 画像データの開始オフセット
    offset = header[10, 4].unpack1("V")
    
    # 幅と高さ (符号付き32ビットリトルエンディアン)
    width = header[18, 4].unpack1("l<")
    height = header[22, 4].unpack1("l<")
    
    # 高さがマイナスの場合はトップダウン(上から下)形式
    is_top_down = height < 0
    height = height.abs

    # 色深度 (1ピクセルあたりのビット数)
    bpp = header[28, 2].unpack1("v")

    unless [1, 24, 32].include?(bpp)
      abort "エラー: 未対応の色深度です(#{bpp}ビット)。1ビットまたは24/32ビットの画像を用意してください。"
    end

    f.seek(offset)
    bitmap_data = []

    if bpp == 1
      # 1ビットBMPの場合: 行のバイト数は4バイト(32ビット)の倍数にパディングされる
      row_bytes = ((width + 31) / 32) * 4
      
      height.times do
        row_data = f.read(row_bytes)
        row_bits = []
        width.times do |x|
          byte_idx = x / 8
          bit_idx = 7 - (x % 8) # MSBから順番にピクセルデータが入る
          bit = (row_data.getbyte(byte_idx) >> bit_idx) & 1 ^ 1 # 0を1、1を0に変換（白を0、黒を1とする場合）
          row_bits << bit
        end
        bitmap_data << row_bits
      end
    else
      # 24/32ビットBMPの場合: RGBの輝度で2値化 (128を閾値とする)
      bytes_per_pixel = bpp / 8
      row_bytes = ((width * bytes_per_pixel + 3) / 4) * 4

      height.times do
        row_data = f.read(row_bytes)
        row_bits = []
        width.times do |x|
          idx = x * bytes_per_pixel
          b = row_data.getbyte(idx)
          g = row_data.getbyte(idx + 1)
          r = row_data.getbyte(idx + 2)
          brightness = (r + g + b) / 3
          # 明るいピクセル(白)を0、暗いピクセル(黒)を1とする場合は条件を調整してください
          row_bits << (brightness > 128 ? 1 : 0)
        end
        bitmap_data << row_bits
      end
    end

    # BMPは通常ボトムアップ（下から上）で保存されるため、配列を反転させる
    bitmap_data.reverse! unless is_top_down

    # 配列名をファイル名から自動生成
    array_name = File.basename(file_path, ".*").gsub(/[^a-zA-Z0-9_]/, '_')
    array_name = "image" if array_name.empty? || array_name =~ /^[0-9]/

    # C言語のコードとして標準出力
    puts "const int #{array_name}_width = #{width};"
    puts "const int #{array_name}_height = #{height};"
    puts "const uint8_t #{array_name}_data[#{height}][#{width}] = {"
    
    bitmap_data.each_with_index do |row, i|
      line = "    {" + row.join(", ") + "}"
      line += "," if i < height - 1
      puts line
    end
    puts "};"
  end
end

if ARGV.empty?
  puts "使用方法: ruby bmp2c.rb <input_file.bmp>"
  exit 1
end

convert_bmp_to_c_array(ARGV[0])