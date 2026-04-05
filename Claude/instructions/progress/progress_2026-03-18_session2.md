# Progress Report: 2026-03-18 Session 2

## Session Summary

ppmax バイナリサーチを実施。156カーネル → 13カーネルまで絞り込み済み。

---

## ppmax バイナリサーチ結果

リファレンス: `log.solver_gpudisabled.ref.txt` (CPU baseline)
比較指標: ppmax at step 60 (user_short.conf, etime=300)
CPU ref ppmax@step60 = **1142.5758**

| Round | テスト | GPU有効カーネル | ppmax@step60 | diff | 結論 |
|-------|--------|----------------|-------------|------|------|
| 1 | PP1 | 003-098 (53個) | 1137.1613 | -0.47% | ← 原因 |
| 1 | PP2 | 102-182 (24個) | 1142.5758 | ≈0 | 無関係 |
| 2 | PP1a | 003-052 (28個) | 1142.5757 | ≈0 | 無関係 |
| 2 | PP1b | 053-098 (25個) | 1137.1613 | -0.47% | ← 原因 |
| 3 | PP1b1 | 053-076 (12個) | 1142.5758 | ≈0 | 無関係 |
| 3 | PP1b2 | 078-098 (13個) | *未テスト* | — | ← 原因（消去法） |

**現在の絞り込み: Kernel 078-098 (13個)**
- 078, 080, 082, 084, 086, 088, 090, 092, 094, 095, 096, 097, 098

### Round 4 準備状況

次の分割:
- **PP1b2a**: 078, 080, 082, 084, 086, 088 (6個) — **ビルド完了、テスト未実行**
- **PP1b2b**: 090, 092, 094, 095, 096, 097, 098 (7個) — compile.conf未作成

---

## 次回セッションでやること

1. **PP1b2a テスト実行**
   - バイナリ: `solver_bsearch_PP1b2a.exe` (ビルド済み)
   - コマンド: `cd test_real && rm -f result/*dmp* result/*mon* result/*geo* && mpirun -np 1 ../solver_bsearch_PP1b2a.exe < user_short.conf > log.bsearch_PP1b2a.txt 2> log.bsearch_PP1b2a.err`
   - compile.conf: 現在の `compile.conf` がPP1b2a用

2. **PP1b2aの結果に応じて**:
   - ppmax≈baseline → 原因はPP1b2b (090-098)
   - ppmax=-0.47% → 原因はPP1b2a (078-088)
   - さらに2分割して絞り込み → 個別カーネル特定 → GPUコード調査・修正

3. **最終検証**: 全156カーネル有効で360stepテスト

---

## ビルド済みバイナリ（本セッション追加分）

| ファイル | 内容 |
|----------|------|
| `solver_bsearch_PP1.exe` | PP1 (003-098 GPU) |
| `solver_bsearch_PP2.exe` | PP2 (102-182 GPU) |
| `solver_bsearch_PP1a.exe` | PP1a (003-052 GPU) |
| `solver_bsearch_PP1b.exe` | PP1b (053-098 GPU) |
| `solver_bsearch_PP1b1.exe` | PP1b1 (053-076 GPU) |
| `solver_bsearch_PP1b2a.exe` | PP1b2a (078-088 GPU) ← テスト未実行 |

## compile.conf保存（本セッション追加分）

| ファイル | 内容 |
|----------|------|
| `compile.conf_bsearch_PP1` | PP1用 |
| `compile.conf_bsearch_PP2` | PP2用 |
| `compile.conf_bsearch_PP1a` | PP1a用 |
| `compile.conf_bsearch_PP1b` | PP1b用 |
| `compile.conf` | **現在PP1b2a用** |

## テストログ（本セッション）

| ファイル | ppmax@step60 | tkemax@step60 | 結論 |
|----------|-------------|--------------|------|
| `log.bsearch_PP1.txt` | 1137.1613 | 3.6295931 | ppmax原因 |
| `log.bsearch_PP2.txt` | 1142.5758 | 3.6295977 | 無関係 |
| `log.bsearch_PP1a.txt` | 1142.5757 | 3.6295991 | 無関係 |
| `log.bsearch_PP1b.txt` | 1137.1613 | 3.6295984 | ppmax原因 |
| `log.bsearch_PP1b1.txt` | 1142.5758 | 3.6295958 | 無関係 |

## 全体の進捗まとめ

- **tkemax問題**: 解決済み（Kernel 204 = more0q.f90、前セッションで修正）
- **ppmax問題**: 078-098 (13カーネル) に絞り込み、次ラウンドで6-7個に
- 残り推定: あと2-3ラウンド（各ラウンド=ビルド+テスト約20分）で個別カーネル特定
