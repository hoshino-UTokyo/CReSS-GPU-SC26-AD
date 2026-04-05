# Progress Report: 2026-03-11 Session 2

## Phase 5: Simulation Integration - TKE Bug Fix

### 成果: swp2nxt.f90 (GPU_321) のTKEバグ修正 — 完了

#### バイナリサーチ結果

| テスト | 無効化ID | TKE結果 |
|--------|----------|---------|
| tke_A1 | 312,313,315,318,320 (前半) | TKE_FAIL |
| tke_A2 | 321,323,327,331,332 (後半) | TKE_PASS |
| tke_B1 | 321,323 | TKE_PASS |
| tke_C1 | 321のみ | TKE_PASS |
| tke_C2 | 323のみ | TKE_FAIL |

→ **GPU_321 (swp2nxt.f90) が唯一の原因**

#### バグの内容

`swp2nxt.f90` の GPU advopt.le.3 ブランチで、**系統的なコピーペーストバグ**。

CPU版では主変数 (u,v,w,pp,ptp,qv等) は3-way swap (`past=current; current=future`)、
しかし aerosol, tracer, TKE, tund は `past=future` のみ（現在値は保持）。

GPU版は全変数に3-way swapを適用していた → **4箇所を修正**:

| 変数 | GPU版(修正前) | CPU版(正) | GPU版(修正後) |
|------|--------------|-----------|--------------|
| aerosol (qaslp) | qaslp=qasl; qasl=qaslf | qaslp=qaslf | qaslp=qaslf |
| tracer (qtp) | qtp=qt; qt=qtf | qtp=qtf | qtp=qtf |
| **TKE (tkep)** | **tkep=tke; tke=tkef** | **tkep=tkef** | **tkep=tkef** |
| tund (tundp) | tundp=tund; tund=tundf | tundp=tundf | tundp=tundf |

#### 検証結果（短縮版 300秒 = 60ステップ）

| 指標 | 修正前 (全156有効) | 修正後 (全156有効) | 基準 |
|------|-------------------|-------------------|------|
| TKE step2 | 1.86e-3 (異常) | 0.0 (正常) | 0.0 |
| tkemax max_reldiff | ~150% | **2.4%** | vs gpudisabled |
| ppmax max_reldiff | 0.57% | 0.57% | vs gpudisabled |
| 発散 | なし | なし | — |

### 全バグ修正の総括

| バグ | ファイル | ID | 症状 | 修正日 |
|------|---------|-----|------|--------|
| ss配列欠落 | steptund.f90 | 311 | ppmax発散 | 2026-03-11 session1 |
| スワップパターン | swp2nxt.f90 | 321 | TKE異常値 | 2026-03-11 session2 |

### tkemax 2.4%差の評価

修正後も tkemax の最大相対差が **2.4%** あるが、既知のFPU精度差で説明可能と判断。

カーネルベンチマーク(Phase 4)で浮動小数点精度の違いによりFAILした5カーネル:

| ID | カーネル | 個別ベンチ誤差 | TKEへの影響 |
|----|---------|--------------|------------|
| 087 | disptke | 66.7% | **TKE散逸を直接計算** (exp/log/sqrt桁落ち) |
| 320 | swadjst | 133% | 飽和調整→熱力学状態→TKEに間接影響 |
| 041 | bruntv | 820% | Brunt-Väisälä振動数→乱流に影響 |
| 296 | siadjst | 0.73% | 飽和判定1ULP差 |
| 056 | cloudcov | 1.9e12% | 雲量(0.0 vs 1.9e-8、実質影響なし) |

**定性的にはFPU精度差で説明可能**だが、定量的な根拠は不十分:
- 上記の誤差率（66.7%等）は配列全体の平均ではなく、**特定の格子点での最大相対誤差**
- 大半の格子点ではCPU/GPU差はもっと小さく、disptkeの66.7%はexp/log/sqrtの
  桁落ちが起きる特定条件の格子点に限られる
- 「disptke 1要素66.7% → 60step蓄積で全体2.4%」と直結させるのは厳密ではない
- ただし、複数カーネル(disptke/swadjst/bruntv等)の精度差が非線形に相互作用して
  蓄積することは確かで、ppmax 0.57%よりTKEの差が大きいのはTKE計算の
  非線形性(exp/log/sqrt)の高さと整合する

**結論: 定性的にはFPU差で説明可能だが、確証にはベースライン比較が必要。**

### Next Steps (次回セッション)

1. **tkemax 2.4%差のベースライン確認**:
   - noacc vs gpudisabled の tkemax差を60ステップ目まで比較
   - gpudisabled ref は361ステップ(1800秒)実行済みなので既存ログで確認可能
   - ベースラインと同等なら許容、超えていれば追加バグの可能性 → 再バイナリサーチ
2. **フル長テスト実行**: user.conf（1800秒）で最終検証
3. **コードクリーンアップ**: テスト用ファイル整理

### ファイル状態

| File | 状態 |
|------|------|
| Src/steptund.f90 | 修正済み（ss配列追加, session1） |
| Src/swp2nxt.f90 | **修正済み（スワップパターン修正, session2）** |
| compile.conf | 全GPU有効 |
| test_real/log.fix321.txt | 全156有効修正後テスト結果（60ステップ） |
| test_real/log.tke_*.txt | バイナリサーチ各ステップのログ |
