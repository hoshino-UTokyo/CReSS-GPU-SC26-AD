# Progress Report: 2026-03-18 Session 5

## Session Summary

GPU (OpenACC) vs CPU (OpenMP 72スレッド) の性能プロファイル比較を実施。

---

## プロファイル条件

- **GPU**: `NVCOMPILER_ACC_TIME=1` による全カーネル計測（360ステップ、全156 GPU有効）
- **CPU**: `Kernel_benchmark/benchmark_results_20260226_215011.log`（72 OpenMPスレッド、per-call平均）
- **ログファイル**: `test_real/log.solver_gpu.20260318151454.err`
- **GPU managed memory** (`-gpu=managed`) 使用のためデータ転送は自動管理

## 全体結果

| 指標 | 値 |
|------|-----|
| GPU総カーネル時間 | 487.6 sec (144ルーチン) |
| CPU推定総時間（マッチ分） | 6,170 sec |
| GPU時間（マッチ分） | 454 sec |
| **全体スピードアップ** | **約13.6倍** |

## GPU時間 Top 15（全体の80%）

| Rank | Routine | CPU/call (ms) | GPU/call (ms) | Speedup | GPU占有率 |
|------|---------|--------------|--------------|---------|----------|
| 1 | pgrad | 45.9 | 3.7 | 12.3x | 11.0% |
| 2 | stepwi | 59.3 | 3.7 | 15.9x | 11.0% |
| 3 | gaussel | 21.1 | 2.7 | 7.7x | 8.1% |
| 4 | diver3d | 34.3 | 2.6 | 13.1x | 7.7% |
| 5 | stepuv | 21.7 | 2.0 | 10.7x | 6.0% |
| 6 | diver2d | 24.3 | 1.8 | 13.3x | 5.4% |
| 7 | advs | 82.4 | 6.5 | 12.6x | 5.3% |
| 8 | advbspi | 13.7 | 1.5 | 9.3x | 4.3% |
| 9 | buoywsi | 20.7 | 1.4 | 14.9x | 4.1% |
| 10 | pgradiv | 17.1 | 1.3 | 13.3x | 3.8% |
| 11 | diverpiv | 8.0 | 0.6 | 12.6x | 3.8% |
| 12 | smoo4s | 52.1 | 4.3 | 12.1x | 3.2% |
| 13 | steppi | 9.0 | 0.9 | 9.5x | 2.8% |
| 14 | diverpih | 7.6 | 0.6 | 12.5x | 1.8% |
| 15 | phy2cnt | 6.4 | 0.6 | 11.2x | 1.8% |

## スピードアップ分布の特徴

### 高速化が大きいカーネル（15x以上）
- `outmxn` 187x — I/O関連、CPUではボトルネックだがGPUでは軽い
- `adjstni` 22x, `setblk` 22x — メモリアクセスパターンがGPU向き
- `swp2nxt` 19x, `stepwi` 16x, `timeflt` 15x

### 典型的なカーネル（10〜15x）
- 大半のカーネルがこの範囲: pgrad(12x), diver3d(13x), diver2d(13x), advs(13x), turbs(14x)
- 3次元ステンシル演算がGPUの大規模並列化に適合

### 高速化が低いカーネル（10x未満）
- `gaussel` 7.7x — Gauss消去法、GPU並列化が本質的に難しい
- `advp` 7.4x, `curveuvw` 6.4x — データアクセスパターンの問題の可能性
- `steps` 2.3x, `depsit` 1.2x, `collect` 4.0x — 計算量が少なくオーバーヘッドが支配的

### GPUが遅いケース
- `setcst3d` **0.47x**（GPUの方が2倍遅い） — 初期化系、呼び出し108回のみ、カーネル起動オーバーヘッドが支配的

## カテゴリ別GPU時間占有率

| カテゴリ | 主要ルーチン | GPU時間合計 | 占有率 |
|---------|------------|-----------|--------|
| 圧力ソルバー | pgrad, stepwi, gaussel, diver3d/2d, steppi, pgradiv, diverpiv/ih | ~270 sec | ~55% |
| 移流系 | advs, advbspi, advuvw, advp | ~57 sec | ~12% |
| 平滑化 | smoo4s, smoo4uvw, smoo4qv | ~22 sec | ~4.5% |
| 浮力 | buoywsi, buoywb | ~21 sec | ~4.3% |
| 乱流 | turbs, turbuvw, turbtke | ~10 sec | ~2% |
| その他 | 残り129ルーチン | ~108 sec | ~22% |

## 最適化の方向性

1. **gaussel (7.7x, GPU時間8.1%)** — 三重対角行列ソルバー。GPU並列化のアルゴリズム改善が最大の伸びしろ（Thomas法 → PCR/CR法など）
2. **setcst3d (0.47x)** — GPU負け。CPU実行にフォールバックするか、カーネル統合で起動オーバーヘッドを削減
3. **advp, curveuvw, steps (6〜7x)** — データアクセスパターンの最適化で改善可能性あり
4. **圧力ソルバー全体（55%）** — カーネル融合やメモリアクセス最適化で大きなインパクト

## Next Steps

1. Nsight Systems / Nsight Compute による詳細プロファイリング
2. gaussel のアルゴリズム検討（PCR法等）
3. managed memory → explicit data management への移行検討
