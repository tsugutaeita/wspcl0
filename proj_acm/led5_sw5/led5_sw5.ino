#define NUM_LED 5
#define LED_00 3
#define LED_01 5
#define LED_02 6
#define LED_03 9
#define LED_04 10
#define NUM_SW 5
#define SW_00 A0
#define SW_01 A1
#define SW_02 A2
#define SW_03 A3
#define SW_04 A4
#define BLINK_INT 100

int led[] = {LED_00, LED_01, LED_02, LED_03, LED_04};
int sw[]  = {SW_00, SW_01, SW_02, SW_03, SW_04};

// 状態を3つ定義
enum LedState { STATE_OFF, STATE_BLINK, STATE_ON };

LedState ledStates[NUM_LED];       // 各LEDの状態
bool lastSwState[NUM_SW];          // スイッチの前回の状態（エッジ検出用）
bool blinkOutput[NUM_LED];         // 点滅時のHIGH/LOW状態
unsigned long lastBlinkMillis = 0; // 点滅用の共通タイマー

void setup() {
  for(int i = 0; i < NUM_LED; i++) {
    pinMode(led[i], OUTPUT);
    digitalWrite(led[i], HIGH); // 初期値：消灯(HIGH)
    ledStates[i] = STATE_OFF;
    blinkOutput[i] = HIGH;
  }
  
  for(int i = 0; i < NUM_SW; i++) {
    // 外付けプルアップ抵抗を使用するため INPUT に設定
    pinMode(sw[i], INPUT);
    lastSwState[i] = HIGH; // 初期値：離されている(HIGH)
  }
}

void loop() {
  unsigned long currentMillis = millis();

  // --- 1. スイッチの入力処理 ---
  for (int i = 0; i < NUM_SW; i++) {
    bool currentSwState = digitalRead(sw[i]);

    // チャタリング防止：状態が変わった瞬間だけ50ms待ち、安定した値を再取得
    if (currentSwState != lastSwState[i]) {
      delay(50);
      currentSwState = digitalRead(sw[i]);
    }

    // 「前回離されていて、今押された瞬間」だけを判定
    if (currentSwState == LOW && lastSwState[i] == HIGH) {
      
      if (ledStates[i] == STATE_BLINK) {
        // 【点滅中だった場合】全消灯して、自分だけ常時点灯
        for (int j = 0; j < NUM_LED; j++) {
          ledStates[j] = STATE_OFF;
          digitalWrite(led[j], HIGH); 
        }
        ledStates[i] = STATE_ON;
        digitalWrite(led[i], LOW); 
        
      } else if (ledStates[i] == STATE_OFF) { 
        // 【修正】消灯中だった場合のみ点滅状態にする (常時点灯中の STATE_ON は無視される)
        
        // もし他に「点滅中」のLEDがあれば、それだけを消灯する（常時点灯はスルー）
        for (int j = 0; j < NUM_LED; j++) {
          if (i != j && ledStates[j] == STATE_BLINK) {
            ledStates[j] = STATE_OFF;
            digitalWrite(led[j], HIGH);
          }
        }
        
        // 自分自身を点滅状態に設定
        ledStates[i] = STATE_BLINK;
        blinkOutput[i] = LOW; // 点灯からスタート
        digitalWrite(led[i], blinkOutput[i]);
      }
    }
    
    // 次回の判定のために状態を記憶
    lastSwState[i] = currentSwState;
  }

  // --- 2. 点滅処理 ---
  if (currentMillis - lastBlinkMillis >= BLINK_INT) {
    lastBlinkMillis = currentMillis;
    
    for (int i = 0; i < NUM_LED; i++) {
      if (ledStates[i] == STATE_BLINK) {
        blinkOutput[i] = !blinkOutput[i]; // HIGH/LOWを反転
        digitalWrite(led[i], blinkOutput[i]);
      }
    }
  }
}
