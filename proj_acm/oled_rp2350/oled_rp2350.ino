#include <Arduino.h>
#include <SPI.h>

// 作成済みのヘッダをインクルード
#include "m_h64w69.h"
#include "shinonome16.h"   // 16x16フォント用
#include "shinonome32.h" // 32x32フォント用

// --- ピン定義 (XIAO RP2350) ---
#define OLED_MOSI D10
#define OLED_SCK  D8
#define OLED_CS   D0
#define OLED_DC   D1
#define OLED_RST  D3

// --- フレームバッファ ---
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

    // RP2350 (arduino-picoコア) 用のSPIピンアサイン
    SPI.setSCK(OLED_SCK);
    SPI.setTX(OLED_MOSI);
    SPI.begin();
    
    SPI.beginTransaction(SPISettings(10000000, MSBFIRST, SPI_MODE0)); // 10MHz

    // SSD1322 基本的な初期化シーケンス
    sendCommand(0xFD); sendData(0x12); // Unlock OLED driver IC
    sendCommand(0xAE);                 // Display OFF

    // カラムアドレスを中央の256ピクセル分(28~91)に設定
    sendCommand(0x15); sendData(0x1C); sendData(0x5B); 
    
    sendCommand(0x75); sendData(0x00); sendData(0x3F); // Set Row Address (0~63)
    sendCommand(0x5C);                 // Enable MCU Write
    sendCommand(0xA0); sendData(0x14); sendData(0x11); // Set Re-map / Dual COM Line Mode
    sendCommand(0xB3); sendData(0xF0); // Front Clock Divider
    sendCommand(0xCA); sendData(0x3F); // MUX Ratio

    // ---------------------------------------------------
    // 【対策3】輝度ムラ（クロストーク）を軽減するための電圧設定
    // ---------------------------------------------------
    // マスターカレント（全体の電流量上限）を少し下げる (0x00~0x0F, デフォルト0x0F)
    sendCommand(0xC7); sendData(0x0F); 
    
    // プレチャージ電圧を引き上げる (0x00~0x1F, デフォルト0x17)
    // ※最大値の 0x1F に設定してムラが改善するか確認します
    sendCommand(0xBB); sendData(0x17); 

    // VCOMH電圧を引き上げる (0x00~0x07, デフォルト0x04)
    // ※最大値の 0x07 に設定してムラが改善するか確認します
    sendCommand(0xBE); sendData(0x04); 
    // ---------------------------------------------------

    sendCommand(0xAF);                 // Display ON
}

// バッファを全クリア (黒で塗りつぶし)
void clearBuffer() {
    memset(frameBuffer, 0, sizeof(frameBuffer));
}

// ---------------------------------------------------------
// ピクセル数ベースの輝度補正テーブル (16段階)
// インデックス0が「点灯数0〜16」、インデックス15が「点灯数241〜256」
// 値(0〜15)が、その段階で適用する輝度です。
// ※ ムラをなくすためのチューニング用配列です。
// ---------------------------------------------------------
const uint8_t brightnessTable[16] = {
    9, // 0:   0〜16粒点灯   (負荷極小 -> あえて輝度を絞る)
    10, // 1:  17〜32粒点灯
    11, // 2:  33〜48粒点灯
    12, // 3:  49〜64粒点灯
    13, // 4:  65〜80粒点灯
    14, // 5:  81〜96粒点灯
    15, // 6:  97〜112粒点灯
    15, // 7: 113〜128粒点灯  (半分点灯)
    15,  // 8: 129〜144粒点灯
    15,  // 9: 145〜160粒点灯
    15,  // 10: 161〜176粒点灯
    15,  // 11: 177〜192粒点灯
    15,  // 12: 193〜208粒点灯
    15,  // 13: 209〜224粒点灯
    15,  // 14: 225〜240粒点灯
    15   // 15: 241〜256粒点灯 (負荷MAX -> MAX輝度)
};

// ---------------------------------------------------------
// ピクセルカウントによる輝度補正関数
// ---------------------------------------------------------
void applyPixelCountCompensation() {
    for (int y = 0; y < DISP_HEIGHT; y++) {
        int litPixelCount = 0;

        // 1. その行の「点灯している粒数」をカウントする
        for (int byteIdx = 0; byteIdx < DISP_WIDTH / 2; byteIdx++) {
            uint8_t data = frameBuffer[y][byteIdx];
            if (data == 0) continue; // 両方黒ならスキップ

            uint8_t left = data >> 4;
            uint8_t right = data & 0x0F;

            // 輝度が0より大きければ「点灯している」とみなす
            if (left > 0) litPixelCount++;
            if (right > 0) litPixelCount++;
        }

        // 2. 点灯数が0なら書き換え不要なので次の行へ
        if (litPixelCount == 0) continue;

        // 3. テーブルのインデックス(0〜15)を算出
        // 例: 1粒なら0、17粒なら1、256粒なら15
        int tableIndex = (litPixelCount - 1) / 16;
        
        // 配列外アクセス防止（念のため）
        if (tableIndex < 0) tableIndex = 0;
        if (tableIndex > 15) tableIndex = 15;

        // 4. テーブルから目標輝度を取得
        uint8_t targetBrightness = brightnessTable[tableIndex];

        // 5. バッファを書き換えて輝度を統一する
        for (int byteIdx = 0; byteIdx < DISP_WIDTH / 2; byteIdx++) {
            uint8_t data = frameBuffer[y][byteIdx];
            if (data == 0) continue;

            uint8_t left = data >> 4;
            uint8_t right = data & 0x0F;

            // 点灯しているピクセルの輝度を、テーブルで引いた値で上書きする
            if (left > 0) left = targetBrightness;
            if (right > 0) right = targetBrightness;

            frameBuffer[y][byteIdx] = (left << 4) | right;
        }
    }
}

// ---------------------------------------------------------
// 真のクロストーク確認用テストパターン (完全版)
// 0〜15のすべての配列インデックスを正確に4行ずつトリガーします
// ---------------------------------------------------------
void generateTestPattern() {
    clearBuffer();

    for (int y = 0; y < DISP_HEIGHT; y++) {
        // 1. 左端に「基準となる縦帯」を描画 (幅4ピクセル: x=0, 1, 2, 3)
        for (int x = 0; x < 4; x++) {
            drawPixel(x, y, 15);
        }

        // 2. 現在の行が16段階のどのステップか (0 〜 15)
        // DISP_HEIGHT=64 なので、4行ごとに1ステップ進む
        int step = y / 4; 

        // 3. このステップで目標とする「合計点灯ピクセル数」
        // 各インデックスの上限値(16, 32, 48 ... 256)を狙う。
        int targetTotalPixels = (step + 1) * 16;
        
        // 【重要】画面幅256pxのうち、隙間(x=4,5,6,7)を空けるため、物理的な最大点灯数は252pxになります。
        // 252pxでもインデックス15（241〜256粒）の範囲内にしっかり収まるため問題ありません。
        if (targetTotalPixels > 252) {
            targetTotalPixels = 252; 
        }

        // 4. 階段部分で点灯すべきピクセル数 = 目標合計 - 縦帯(4)
        int staircasePixels = targetTotalPixels - 4;

        // 5. X=8 から階段を描画
        for (int x = 8; x < 8 + staircasePixels; x++) {
            if (x < DISP_WIDTH) { // バッファあふれ防止
                drawPixel(x, y, 15);
            }
        }
    }
}

// バッファの内容をOLEDに転送
void updateDisplay() {
    sendCommand(0x15); sendData(0x1C); sendData(0x5B);
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

    // 偶数Xは上位4bit、奇数Xは下位4bit
    if (x % 2 == 0) {
        frameBuffer[y][byteIdx] = (frameBuffer[y][byteIdx] & 0x0F) | (color << 4);
    } else {
        frameBuffer[y][byteIdx] = (frameBuffer[y][byteIdx] & 0xF0) | color;
    }
}

// 画像の描画 (左端に描画)
void drawImage() {
    for (int y = 0; y < m_mono_H64W69_height; y++) {
        for (int x = 0; x < m_mono_H64W69_width; x++) {
            if (m_mono_H64W69_data[y][x] == 1) {
                drawPixel(x, y, 10); // ★ムラ対策で輝度10推奨
            }
        }
    }
}

// ====================================================
// 16x16フォント描画関数群
// ====================================================
void drawChar16(int x, int y, uint32_t utf8_code) {
    for (int i = 0; i < shinonome16_num_glyphs; i++) {
        if (shinonome16_font[i].utf8_code == utf8_code) {
            for (int row = 0; row < 16; row++) {
                uint8_t b1 = shinonome16_font[i].bitmap[row * 2];
                uint8_t b2 = shinonome16_font[i].bitmap[row * 2 + 1];

                for (int col = 0; col < 8; col++) {
                    if ((b1 >> (7 - col)) & 1) drawPixel(x + col, y + row, 15);
                    if ((b2 >> (7 - col)) & 1) drawPixel(x + 8 + col, y + row, 15);
                }
            }
            break;
        }
    }
}

void drawString16(int x, int y, const char* str) {
    int cursor_x = x;
    int cursor_y = y;
    int i = 0;

    while (str[i] != '\0') {
        uint8_t c = str[i];
        uint32_t code = 0;
        int char_width = 16; 

        if ((c & 0x80) == 0x00) {
            code = c;
            i += 1;
            char_width = 8; // 半角は8px進む
        } else if ((c & 0xE0) == 0xC0) {
            code = (c << 8) | (uint8_t)str[i+1];
            i += 2;
        } else if ((c & 0xF0) == 0xE0) {
            code = ((uint32_t)c << 16) | ((uint32_t)(uint8_t)str[i+1] << 8) | (uint8_t)str[i+2];
            i += 3;
        } else {
            i++;
            continue;
        }

        if (code == '\n') {
            cursor_x = x;
            cursor_y += 16;
            continue;
        }

        drawChar16(cursor_x, cursor_y, code);
        cursor_x += char_width; 
        
        if (cursor_x + char_width > DISP_WIDTH) {
            cursor_x = x;
            cursor_y += 16;
        }
    }
}

// ====================================================
// 32x32フォント描画関数群
// ====================================================
void drawChar32(int x, int y, uint32_t utf8_code) {
    for (int i = 0; i < shinonome32_num_glyphs; i++) {
        if (shinonome32_font[i].utf8_code == utf8_code) {
            for (int row = 0; row < 32; row++) {
                uint8_t b1 = shinonome32_font[i].bitmap[row * 4];     // 1~8px
                uint8_t b2 = shinonome32_font[i].bitmap[row * 4 + 1]; // 9~16px
                uint8_t b3 = shinonome32_font[i].bitmap[row * 4 + 2]; // 17~24px
                uint8_t b4 = shinonome32_font[i].bitmap[row * 4 + 3]; // 25~32px

                for (int col = 0; col < 8; col++) {
                    if ((b1 >> (7 - col)) & 1) drawPixel(x + col,      y + row, 15);
                    if ((b2 >> (7 - col)) & 1) drawPixel(x + 8 + col,  y + row, 15);
                    if ((b3 >> (7 - col)) & 1) drawPixel(x + 16 + col, y + row, 15);
                    if ((b4 >> (7 - col)) & 1) drawPixel(x + 24 + col, y + row, 15);
                }
            }
            break;
        }
    }
}

void drawString32(int x, int y, const char* str) {
    int cursor_x = x;
    int cursor_y = y;
    int i = 0;

    while (str[i] != '\0') {
        uint8_t c = str[i];
        uint32_t code = 0;
        int char_width = 32; // 全角デフォルト

        if ((c & 0x80) == 0x00) {
            code = c;
            i += 1;
            char_width = 16; // 半角は16px進む
        } else if ((c & 0xE0) == 0xC0) {
            code = (c << 8) | (uint8_t)str[i+1];
            i += 2;
        } else if ((c & 0xF0) == 0xE0) {
            code = ((uint32_t)c << 16) | ((uint32_t)(uint8_t)str[i+1] << 8) | (uint8_t)str[i+2];
            i += 3;
        } else {
            i++;
            continue;
        }

        if (code == '\n') {
            cursor_x = x;
            cursor_y += 32;
            continue;
        }

        drawChar32(cursor_x, cursor_y, code);
        cursor_x += char_width; 
        
        if (cursor_x + char_width > DISP_WIDTH) {
            cursor_x = x;
            cursor_y += 32;
        }
    }
}

// ====================================================
// メイン処理
// ====================================================
void setup() {
    Serial.begin(115200);
    
    initSSD1322();
    clearBuffer();
    drawImage();
    
    // 【描画テスト】
    // Y=0〜31 の領域に 32x32 フォントを描画
    // drawString32(72, 0, "　会議中　");
    // drawString32(72, 32, " 左ドアへ ");

    // Y=32〜 の領域に 16x16 フォントを描画
    // drawString16(72, 32, "RP2350 16x16");
    // drawString16(72, 48, "SSD1322 OLED");
    // drawString16(72,  0, "AUTHENTICATION SUCCESS");
    // drawString16(72, 16, "ACTIVATION SUCCESS");
    // drawString16(72, 32, "UNLOCK COMPLETE");
    // drawString16(72, 48, "OMEGA PSYCOMMU ACTIVE");
    drawString32(72, 0, "直流が命。");
    drawString32(80, 32, "   和田和千");
    
    // ★ テストパターンを描画バッファに生成
    // generateTestPattern();
    
    // ★ 描画が終わったバッファに対して、ピクセル数ベースの補正をかける
    applyPixelCountCompensation();

    updateDisplay();
}

void loop() {
}