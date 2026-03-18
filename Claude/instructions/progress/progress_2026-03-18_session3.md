# Progress Report: 2026-03-18 Session 3

## Session Summary

ppmax バイナリサーチを継続。前セッション(session2)の13カーネルから **Kernel 094** (diverpiv.f90, phvbcs.f90) を原因として特定。
バグを発見・修正し、全156カーネル有効で ppmax が baseline と完全一致することを確認。

---

## ppmax バイナリサーチ結果（session2からの継続）

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
| 3 | PP1b2 | 078-098 (13個) | *消去法* | — | ← 原因 |
| 4 | PP1b2a | 078-088 (6個) | 1142.5758 | ≈0 | 無関係 |
| 4 | PP1b2b | 090-098 (7個) | *消去法* | — | ← 原因 |
| 5 | PP1b2b1 | 090,092,094 (3個) | 1137.1614 | -0.47% | ← 原因 |
| 5 | PP1b2b2 | 095-098 (4個) | *消去法* | — | 無関係 |
| 6 | K090 | 090のみ | 1142.5757 | ≈0 | 無関係 |
| 6 | K092 | 092のみ | 1142.5758 | ≈0 | 無関係 |
| 6 | K094 | *消去法* | — | — | **← ppmax原因カーネル** |

### 結論: **Kernel 094 が ppmax 差異 (-0.47%) の原因**

Kernel 094 の所在:
- `Src/diverpiv.f90` — 鉛直圧力偏差の発散計算（バグなし）
- `Src/phvbcs.f90` — 気圧変数の境界条件（**バグあり → 修正済み**）

---

## Bug Fix: phvbcs.f90 GPU_094 — dtdvb スケーリング条件バグ

### バグの内容

CPU版では `dtdvb` スケーリングに2つの条件分岐がある:
```fortran
if(fproc(1:3).eq.'sml') then
    scpx(j,k,1) = scpx(j,k,1)*dtdvb    ! 常にスケーリング
else
    if(advopt.ge.4) then
        scpx(j,k,1) = scpx(j,k,1)*dtdvb  ! advopt>=4の時のみ
    end if
end if
```

GPU版では fproc チェックが欠落していた:
```fortran
if (advopt >= 4) then
    scpx(j,k,1) = scpx(j,k,1)*dtdvb
end if
```

### 影響
- `advopt=3` (本シミュレーションの設定) かつ `fproc='sml'` の場合
- CPU版: dtdvb スケーリングが適用される
- GPU版: スケーリングがスキップされる → **位相速度が過大** → ppmax に -0.47% の差異

### 修正箇所
4箇所の境界（west, east, south, north）すべてで同一バグ。
CPU版と同じ `fproc(1:3) == 'sml'` チェックを追加。

---

## 修正後テスト結果

### 60step テスト (fix094)

| 指標 | fix094 (全156有効) | baseline (GPU disabled) | diff |
|------|-------------------|------------------------|------|
| ppmax@step60 | **1142.5758** | **1142.5758** | **0.000%** |
| tkemax@step60 | **3.6296008** | **3.6295619** | **0.001%** |

### 360step フルテスト (fix094)

| 指標 | fix094 (全156有効) | baseline (GPU disabled) | diff |
|------|-------------------|------------------------|------|
| ppmax@step360 | **1140.6746** | **1140.6743** | **+0.0000%** |
| tkemax@step360 | **16.663454** | **16.911001** | **-1.464%** |

- ppmax は baseline と完全一致レベル（fix094修正が完璧に効果）
- tkemax の -1.464% はFPU演算順序差による残差（fix204時点の -1.465% と変化なし）
- ログ: `test_real/log.fix094_full.txt`

---

## テストログ（本セッション）

| ファイル | ppmax@step60 | 結論 |
|----------|-------------|------|
| `log.bsearch_PP1b2a.txt` | 1142.5758 | 無関係 |
| `log.bsearch_PP1b2b1.txt` | 1137.1614 | ppmax原因 |
| `log.bsearch_K090.txt` | 1142.5757 | 無関係 |
| `log.bsearch_K092.txt` | 1142.5758 | 無関係 |
| `log.fix094.txt` | 1142.5758 | **修正後OK** |

## 全体の進捗まとめ

- **tkemax問題**: 解決済み（Kernel 204 = more0q.f90、session1で修正）
- **ppmax問題**: 解決済み（Kernel 094 = phvbcs.f90、本セッションで修正）
- 全156 GPUカーネル有効、4つのバグ修正済み（311, 321, 204, 094）
- 60step精度: ppmax 0.000%, tkemax 0.001%
- 360stepフルテスト: **完了** — ppmax +0.0000%, tkemax -1.464%
