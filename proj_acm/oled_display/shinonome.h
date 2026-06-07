#ifndef SHINONOME_H
#define SHINONOME_H

#include <stdint.h>

typedef struct {
    uint32_t utf8_code;  // UTF-8のバイト列を結合した数値 (例: 'A' -> 0x41, 'あ' -> 0xE38182)
    uint8_t bitmap[32];  // 16x16ピクセル = 256ビット = 32バイト (横8ビット×2バイトが1行、それを16行)
} ShinonomeGlyph;

// 検索用に文字コード順でソートされていると二分探索が使えてベターですが、
// 今回は実装をシンプルにするため線形探索を想定します。
const int shinonome_num_glyphs = 3; 

const ShinonomeGlyph shinonome_font[] = {
    // [注意] バックスラッシュ(0x5C)の定義例
    { 0x0000005C, { 0x80, 0x00, 0x40, 0x00, /* 略... 斜めの線のデータ */ } },
    { 0x00000041, { /* 'A' のデータ */ } },
    { 0x00E38182, { /* 'あ' のデータ */ } }
};

#endif