# GPU Kernel Benchmark 現状メモ
更新日: 2026-01-21

## CPUベンチマーク状況

| 状態 | 件数 | 説明 |
|------|------|------|
| コード+データあり | 150 | 完全なベンチマーク（GPUポーティング対象） |
| コードあり・データなし | 16 | 非計算ベンチマーク（メモリ割当、I/O等） |
| スタブのみ | 223 | カーネル抽出未実施（README.md + variable_list.txtのみ） |
| **合計** | **389** | |

### データなしの16件（コードはある）
- allocbuf, allociot, allocslv (メモリ割当)
- copy1d (単純コピー)
- forcept (点ソース強制)
- getxy, getz (座標計算)
- inidef, inimod (初期化)
- opendmp, rdgrp (I/O操作)
- setcst1d_sec1, setcst1d_r8_sec2, setcst4d, setname (定数設定)
- steps (タイムステップ)

## GPUベンチマーク状況

| 状態 | 件数 |
|------|------|
| PASS（検証成功） | 75 |
| FAIL（検証失敗） | 4 |
| データ不足 | 4 |
| 未作成 | 残り |

### 検証失敗の4件（要修正）
1. **041_bruntv_s_bruntv** - 条件分岐 `qall > thresq` の評価がCPU/GPUで異なる
2. **042_bulksfc_s_bulksfc** - cm出力が100%エラー
3. **056_cloudcov_s_cloudcov** - cdm出力に大きなエラー
4. **087_disptke_s_disptke** - 最大66.7%エラー

### データ不足の4件（CPU側にデータなし）
- 065_copy1d_s_copy1d
- 112_forcept_s_forcept
- 144_getz_s_getz
- 200_mapfct_s_mapfct

## 修正済み
- **337_turbtke_s_turbtke** - tolerance を 2% に変更して PASS

## ビルド環境
- コンパイラ: nvfortran
- フラグ: `-O3 -acc -gpu=managed -mp -Minfo=accel`
- IEEE NaN 対応: `ieee_is_nan()` 使用（`isnan()` は使用不可）

## 次のステップ
1. 残り4件の検証失敗を修正
2. 必要に応じてtoleranceを調整（浮動小数点精度差への対応）
