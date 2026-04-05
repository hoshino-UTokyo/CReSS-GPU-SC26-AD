# Progress Report: 2026-03-06 Session 2

## Phase 5: Simulation Integration - Binary Search (Corrected)

### Summary

前回セッションのバイナリサーチに重大な問題を発見し修正。前回は「orig 19」カーネル（177, 217, 218, 220, 232, 234, 238, 281-283, 286-287, 327, 380-387）が常に無効化された状態でバイナリサーチを行っていたが、Run 10以降でこれらを有効化してしまい結果が矛盾していた。今回は156個全カーネルを対象に正しいバイナリサーチをやり直し、問題を20個のカーネル（行118-137）に絞り込んだ。

### 前回セッションの問題

前回のバイナリサーチ（Run 4-9）では `compile.conf_bsearch_F1` をベースに使用しており、19個の「orig」カーネルが常に無効化されていた。Run 10以降（本セッション初期含む）ではこれらを有効化してテストしたため、全てFAILとなり混乱を招いた。

### 今回のバイナリサーチ結果

`/tmp/all_gpu_ids.txt` に156カーネルIDを行番号順に格納。バイナリサーチスクリプト `bsearch_test.sh` を作成。

| Run | Disabled Lines | Disabled Range | Result | Conclusion |
|-----|---------------|---------------|--------|------------|
| A1 | 1-78 | 003-190 | **FAIL** (step ~55) | Problem in lines 79-156 (200-387) |
| A2 | 79-156 | 200-387 | **PASS** | Confirms problem in 79-156 |
| B1 | 79-117 | 200-291 | **FAIL** (step ~55) | Problem in lines 118-156 (294-387) |
| B2 | 118-156 | 294-387 | **PASS** | Confirms problem in 118-156 |
| C1 | 118-137 | 294-332 | **PASS** | Problem IS in 118-137 |
| C2 | 138-156 | 333-387 | **FAIL** (step ~55) | Confirms problem in 118-137 |

**バイナリサーチロジック**: 「disable X → PASS」は、Xの範囲がGPU上にいなくても問題ない（Xに問題カーネルを含む可能性がある）。「disable X → FAIL」は、X以外の範囲に問題がある。

### 現在の絞り込み: 20カーネル（lines 118-137）

| Line | Kernel ID | Source File |
|------|-----------|-------------|
| 118 | 294 | sheartke.f90 |
| 119 | 295 | shedding.f90 |
| 120 | 296 | siadjst.f90 |
| 121 | 302 | smoo4qv.f90 |
| 122 | 303 | smoo4s.f90 |
| 123 | 304 | smoo4uvw.f90 |
| 124 | 305 | sndwave.f90 |
| 125 | 308 | steppi.f90 |
| 126 | 310 | steps.f90 |
| 127 | 311 | steptund.f90 |
| 128 | 312 | steptund.f90 |
| 129 | 313 | stepuv.f90 |
| 130 | 315 | stepwi.f90 |
| 131 | 318 | strsten.f90 |
| 132 | 320 | swadjst.f90 |
| 133 | 321 | swp2nxt.f90 |
| 134 | 323 | termblk.f90 |
| 135 | 327 | timeflt.f90 |
| 136 | 331 | totalqwi.f90 |
| 137 | 332 | totals.f90 |

### Next Steps (Priority Order)

1. **D1テスト**: lines 118-127 (294-311) を無効化して実行 → PASS/FAILで前半/後半を判定
2. **さらにバイナリサーチ**して、最終的に1-3個のカーネルに絞る
3. **ソースコード修正**: CPU版とGPU版を比較して、ループ範囲・変数スコープ・レースコンディション等を修正
4. **修正後、全156カーネルでフル360ステップ検証**: a_1, a_2 < 1.0e-4 を確認

### Tools/Files Created

| File | Purpose |
|------|---------|
| `bsearch_test.sh` | バイナリサーチ自動化スクリプト（clean+build+run+判定） |
| `/tmp/all_gpu_ids.txt` | 156カーネルIDの完全リスト（行番号でバイナリサーチ） |
| `test_real/log.bsearch_run{10-14}.txt` | 本セッション初期のテストログ（orig-19問題発覚前） |
| `test_real/log.bsearch_{A1,A2,B1,B2,C1,C2}.txt` | 新バイナリサーチのログ |

### 注意事項

- `/tmp/all_gpu_ids.txt` は一時ファイル。次セッションで再生成が必要:
  ```bash
  grep -rhoP 'DISABLE_GPU_\d+' Src/*.f90 | sort -t_ -k3 -n | uniq > /tmp/all_gpu_ids.txt
  ```
- `bsearch_test.sh` はこのファイルに依存
- compile.confの現在の状態: bsearch_C2（lines 138-156 disabled）のまま
- Interactive jobで直接GPU実行が可能な環境で作業すること
