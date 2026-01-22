# GPU Benchmark FAILED 解析結果

日付: 2026-01-22

## 概要

7つのFAILEDベンチマークを解析した結果、以下の2種類の問題に分類される：

1. **浮動小数点精度の問題** (5件) - GPU移植自体は正しいが、CPU/GPU間のFPU実装差による
2. **データ生成の問題** (2件) - ベンチマークデータ自体に欠陥がある

---

## 1. 浮動小数点精度の問題

### 1.1 296_siadjst_s_siadjst (0.73% error)

**症状**: 境界条件 `qi > dqi` の判定でCPUとGPUが異なる結果

**原因**:
- CPU: qi = X (bit pattern: ...F7)
- GPU: qi = X-1ULP (bit pattern: ...F6)
- 1 ULP (Unit in Last Place) の差で異なる分岐に入る

**詳細**: 飽和調整計算で、qi と dqi がほぼ等しい場合、微小な計算誤差により異なる分岐パスを通る

---

### 1.2 041_bruntv_s_bruntv (820% error)

**症状**: 境界条件 `t <= tlow` (tlow=233.16K) の判定差

**原因**:
- CPU: t = 233.16002 → t <= tlow = False
- GPU: t = 233.16000 → t <= tlow = True
- 1 ULP差で異なる分岐に入り、潜熱項の加算有無が変わる

**詳細**:
```fortran
t3d(i,j,k) = pt(i,j,k) * exp(rddvcp * log(p0iv * (pbr(i,j,k) + pp(i,j,k))))
if (t3d(i,j,k) <= tlow) then
  a(i,j,k) = a(i,j,k) + (lf0 + cwmci * (t3d(i,j,k) - t0))  ! この行の実行有無が変わる
end if
```

---

### 1.3 056_cloudcov_s_cloudcov (1.9e12 error)

**症状**: 極端に大きな相対誤差

**原因**:
- GPU: cdm = 1.9e-8 (微小な非ゼロ値)
- Ref: cdm = 0.0 (正確なゼロ)
- 相対誤差 = 1.9e-8 / 1e-20 = 1.9e12

**詳細**: 境界条件 `zph >= 4800m` の判定で、GPUが微小に異なる結果を計算し、正確にゼロになるべき値が微小な非ゼロ値になる

---

### 1.4 087_disptke_s_disptke (66.7% error)

**症状**: Catastrophic cancellation (壊滅的キャンセル)

**原因**:
- Input: 4.51763335E-04
- CPU term: 4.51763422E-04 → result = -8.73E-11
- GPU term: 4.51763350E-04 → result = -1.46E-10
- term計算の微小差 (~7E-10) がほぼ同じ値の引き算で増幅

**詳細**:
```fortran
ln = (priv - 1.0) * exp(oned3 * log(ds308 * rmf * jcb)) + eps
tkefrc = tkefrc - coefficient * rst * tke * sqrt(tke) / ln
```
exp, log, sqrt の GPU実装差が最終結果に大きく影響

---

### 1.5 320_swadjst_s_swadjst (133% error)

**症状**: qc (雲水混合比) の境界条件判定差

**原因**:
- GPU: qc = 1.93e-9
- Ref: qc = 8.26e-10
- 境界条件 `qc > dqc` の判定で異なる分岐

**詳細**: 飽和調整の2回反復計算で、dqc計算のわずかな差が最終結果に影響

---

## 2. データ生成の問題

### 2.1 042_bulksfc_s_bulksfc (100% error)

**症状**: CPU版も同じエラーでFAIL (1,612,808 errors)

**原因**: params.txt に必要なパラメータが欠損

**欠損パラメータ**:
- wkappa (水面のvon Kármán定数)
- prnumg, prnumw (地表・水面のPrandtl数)
- icz0m, icz0h (氷面の粗度長)
- rmg, rms, rhg, rhs (安定度関数パラメータ)
- oned3 (1/3)
- cc (定数)

**対策**: ダンプデータ再生成が必要

---

### 2.2 239_phycood_s_phycood_sec3 (12500 error)

**症状**: CPU版も同じエラーでFAIL (103,671,720 errors)

**原因**: params.txt のパラメータ値が不正

**問題のパラメータ**:
```
nkm1 = 0    ← 本来 nk-1 = 127 であるべき
nkm2 = 0    ← 本来 nk-2 = 126 であるべき
htuiv = 0.0
htuivz = 0.0
```

**影響**: `zph(i,j,nkm1)` で配列範囲外アクセス (zphは1:nk)

**対策**: ダンプデータ再生成が必要

---

## 根本原因の技術的説明

### CPU/GPU間のFPU実装差

1. **FMA (Fused Multiply-Add)**: GPUはFMA命令を積極的に使用し、中間丸めを省略
2. **超越関数 (exp, log, sqrt)**: 実装アルゴリズムが異なり、最下位ビットが異なる
3. **演算順序**: コンパイラ最適化により演算順序が変わる可能性

### 問題が顕在化するパターン

1. **境界条件の判定**: `if (x > threshold)` で閾値近傍の値
2. **Catastrophic cancellation**: ほぼ等しい値の引き算
3. **ゼロ近傍の相対誤差**: 分母が極小値の場合

---

## 対策オプション

### コンパイラフラグによる対策 (ユーザーが拒否)
- `-Kieee` または `-gpu=nofma`: IEEE準拠モードで精度向上
- 副作用: 性能低下

### 許容誤差の調整 (ユーザーが拒否)
- tolerance を 1e-5 から 1e-3 程度に緩和
- 副作用: 真のバグを見逃す可能性

### コード修正による対策 (限定的)
- 境界条件にイプシロン余裕を追加: `if (x > threshold - eps)`
- 副作用: 物理的意味が変わる可能性

### データ再生成
- 042_bulksfc, 239_phycood_sec3 はダンプデータの再生成が必要

---

## 結論

- 5件の浮動小数点精度問題は、GPU移植自体は正しく、CPU/GPU間の本質的な差異
- 2件のデータ問題は、ベンチマークデータ生成時の問題
- コンパイラフラグや許容誤差を変更せずに解決することは困難
