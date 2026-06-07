#include <Arduino.h>
#include <SPI.h>

// 作成済みのヘッダをインクルード
#include "m_h64w69.h"
#include "shinonome.h"

// --- ピン定義 (XIAO ESP32-C3) ---
#define OLED_MOSI D10
#define OLED_SCK  D8
#define OLED_CS   D7
#define OLED_DC   D6
#define OLED_RST  D5

// --- フレームバッファ ---
// SSD1322は4bitグレースケール。1バイトで2ピクセル(横)を表現。
// 幅256ピクセル = 128バイト, 高さ64ピクセル
#define DISP_WIDTH  256
#define DISP_HEIGHT 64
uint8_t frameBuffer[DISP_HEIGHT][DISP_WIDTH / 2];

// コマンド送信関数
void sendCommand(uint8_t cmd) {
    digitalWrite(OLED_DC, LOW);
    digitalWrite(OLED_CS, LOW);
    SPI.transfer(cmd);
    digitalWrite(OLED_CS, HIGH);
}

// データ送信関数
void sendData(uint8_t data) {
    digitalWrite(OLED_DC, HIGH);
    digitalWrite(OLED_CS, LOW);
    SPI.transfer(data);
    digitalWrite(OLED_CS, HIGH);
}

// SSD1322初期化
void initSSD1322() {
    pinMode(OLED_CS, OUTPUT);
    pinMode(OLED_DC, OUTPUT);
    pinMode(OLED_RST, OUTPUT);

    digitalWrite(OLED_CS, HIGH);

    // ハードウェアリセット
    digitalWrite(OLED_RST, LOW);
    delay(10);
    digitalWrite(OLED_RST, HIGH);
    delay(10);

    SPI.begin(OLED_SCK, -1, OLED_MOSI, OLED_CS);
    SPI.beginTransaction(SPISettings(10000000, MSBFIRST, SPI_MODE0)); // 10MHz

    // SSD1322 基本的な初期化シーケンス
    sendCommand(0xFD); sendData(0x12); // Unlock OLED driver IC
    sendCommand(0xAE);                 // Display OFF
    sendCommand(0x15); sendData(0x00); sendData(0x7F); // Set Column Address (0~127 = 256 pixels)
    sendCommand(0x75); sendData(0x00); sendData(0x3F); // Set Row Address (0~63)
    sendCommand(0x5C);                 // Enable MCU Write
    sendCommand(0xA0); sendData(0x14); sendData(0x11); // Set Re-map / Dual COM Line Mode
    sendCommand(0xB3); sendData(0xF0); // Front Clock Divider
    sendCommand(0xCA); sendData(0x3F); // MUX Ratio
    sendCommand(0xAF);                 // Display ON
}

// バッファを全クリア (黒で塗りつぶし)
void clearBuffer() {
    memset(frameBuffer, 0, sizeof(frameBuffer));
}

// バッファの内容をOLEDに転送
void updateDisplay() {
    sendCommand(0x15); sendData(0x00); sendData(0x7F);
    sendCommand(0x75); sendData(0x00); sendData(0x3F);
    sendCommand(0x5C); // Write RAM command

    digitalWrite(OLED_DC, HIGH);
    digitalWrite(OLED_CS, LOW);
    for (int y = 0; y < DISP_HEIGHT; y++) {
        for (int x = 0; x < (DISP_WIDTH / 2); x++) {
            SPI.transfer(frameBuffer[y][x]);
        }
    }
    digitalWrite(OLED_CS, HIGH);
}

// ピクセル描画関数 (color: 0=黒 ~ 15=白)
void drawPixel(int x, int y, uint8_t color) {
    if (x < 0 || x >= DISP_WIDTH || y < 0 || y >= DISP_HEIGHT) return;
    
    int byteIdx = x / 2;
    color &= 0x0F; // 4bitにクリップ

    // 偶数Xは上位4bit、奇数Xは下位4bit (OLEDのRemap設定による)
    if (x % 2 == 0) {
        frameBuffer[y][byteIdx] = (frameBuffer[y][byteIdx] & 0x0F) | (color << 4);
    } else {
        frameBuffer[y][byteIdx] = (frameBuffer[y][byteIdx] & 0xF0) | color;
    }
}

// Ruby出力画像の描画 (左端に描画)
void drawImage() {
    // ヘッダ名に合わせて m_h64w69_width / m_h64w69_height / m_h64w69_data を使用
    for (int y = 0; y < m_h64w69_height; y++) {
        for (int x = 0; x < m_h64w69_width; x++) {
            // 値が1(描画対象)なら輝度15(真っ白)で描画
            if (m_h64w69_data[y][x] == 1) {
                drawPixel(x, y, 15);
            }
        }
    }
}

// UTF-8コードをキーにして東雲フォント(1文字)を描画
void drawChar(int x, int y, uint32_t utf8_code) {
    // 線形探索でフォントを探す
    for (int i = 0; i < shinonome_num_glyphs; i++) {
        if (shinonome_font[i].utf8_code == utf8_code) {
            // 見つかったら 16x16 ピクセルを描画
            for (int row = 0; row < 16; row++) {
                uint8_t b1 = shinonome_font[i].bitmap[row * 2];     // 左8px
                uint8_t b2 = shinonome_font[i].bitmap[row * 2 + 1]; // 右8px

                for (int col = 0; col < 8; col++) {
                    if ((b1 >> (7 - col)) & 1) drawPixel(x + col, y + row, 15);
                    if ((b2 >> (7 - col)) & 1) drawPixel(x + 8 + col, y + row, 15);
                }
            }
            break;
        }
    }
}

// UTF-8文字列を解析して描画 (16x16フォント前提)
void drawString(int x, int y, const char* str) {
    int cursor_x = x;
    int cursor_y = y;
    int i = 0;

    while (str[i] != '\0') {
        uint8_t c = str[i];
        uint32_t code = 0;
        int char_width = 16; // 全角デフォルト

        if ((c & 0x80) == 0x00) {
            // 1バイト文字 (ASCII)
            code = c;
            i += 1;
            char_width = 8; // 半角の場合は幅を8にする等、必要に応じて調整
        } else if ((c & 0xE0) == 0xC0) {
            // 2バイト文字
            code = (c << 8) | (uint8_t)str[i+1];
            i += 2;
        } else if ((c & 0xF0) == 0xE0) {
            // 3バイト文字 (一般的な日本語)
            code = ((uint32_t)c << 16) | ((uint32_t)(uint8_t)str[i+1] << 8) | (uint8_t)str[i+2];
            i += 3;
        } else {
            // 4バイト以上はスキップ(未対応とする)
            i++;
            continue;
        }

        // 改行処理が必要であればここに実装
        if (code == '\n') {
            cursor_x = x;
            cursor_y += 16;
            continue;
        }

        drawChar(cursor_x, cursor_y, code);
        cursor_x += char_width; // 次の文字へ
        
        // 画面端の折り返し処理
        if (cursor_x + char_width > DISP_WIDTH) {
            cursor_x = x;
            cursor_y += 16;
        }
    }
}

void setup() {
    Serial.begin(115200);
    
    // 1. OLEDの初期化
    initSSD1322();
    
    // 2. バッファのクリア
    clearBuffer();
    
    // 3. 画像の描画 (左端の 69x64 領域)
    drawImage();
    
    // 4. フォントの描画 (画像と被らないようにX=72あたりから開始)
    // エスケープシーケンスに注意: C:\ を表示するには C:\\ と書く
    drawString(72, 8, "ESP32-C3");
    drawString(72, 24, "SSD1322 OLED");
    drawString(72, 40, "東雲ﾌｫﾝﾄ \\0x5C");

    // 5. 画面に反映
    updateDisplay();
}

void loop() {
    // 今回は静的表示のためループ内は空。
    // 追記したい場合は、ここで drawString(...) を呼び、
    // updateDisplay() を実行すればバッファの差分（追記内容）が反映されます。
}