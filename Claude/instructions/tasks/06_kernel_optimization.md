# Phase 6: Kernel Optimization (nsys Roofline Analysis)

[Back to Main](../main.md)

## Objective

nsys (NVIDIA Nsight Systems) によるルーフライン解析に基づき、個別GPUカーネルの性能最適化を行う。最適化は `Kernel_benchmark_gpu_opt/` ディレクトリで独立に実施し、検証後に `Src/` へ統合する。

---

## Prerequisites

- [ ] Phase 5 (Simulation Integration) が完了し、全156カーネルが統合済み
- [ ] 全カーネルの精度検証が合格済み
- [ ] nsys がインストール済み (`nsys --version`)

---

## 最適化の制約条件

### 許可される最適化

| カテゴリ | 例 |
|----------|-----|
| ループ構造の変更 | collapse の追加・削除、ループ分割・融合 |
| OpenACC ディレクティブの調整 | gang/worker/vector の明示指定、tile 節 |
| メモリアクセスパターンの改善 | ループ順序の入れ替え（k-j-i → i-j-k 等）、一時配列の導入 |
| 算術の最適化 | 冗長計算の除去、除算→逆数乗算、条件分岐の単純化 |
| ループ内サブルーチンのインライン化 | `!$acc routine seq` 呼び出しの手動インライン |
| スカラ一時変数の活用 | 配列アクセス回数の削減 |

### 禁止される最適化

以下の最適化は**他カーネルとの依存関係を生じさせるため禁止**:

| カテゴリ | 理由 |
|----------|------|
| `!$acc async` 節による非同期実行 | カーネル間の実行順序依存が発生。async を使う場合は全体のタイムステップ設計を変更する必要がある |
| データ構造の変更 (AoS↔SoA 等) | 他のカーネルが同じ配列を参照しているため、全体の変更が必要になる |
| `!$acc data` / `!$acc enter data` によるデータ配置制御 | Unified Memory 前提の設計と矛盾し、他カーネルのデータ可視性に影響する |
| 配列の次元数・形状・アロケーション変更 | 呼び出し元や他カーネルとのインターフェースが変わる |
| MPI 通信パターンの変更 | ドメイン分割・ハロ交換は全体設計の一部 |
| カーネル間でのデータの事前計算・キャッシュ | 他カーネルの入出力契約が変わる |

**原則**: 各カーネルのサブルーチンの入出力引数（個数・型・形状）は変更しない。内部実装のみ最適化する。

---

## ディレクトリ構成

```
Kernel_benchmark_gpu_opt/
├── Makefile.common               # 共通コンパイル設定（gpu版からコピー）
├── build_all.sh                  # 全ベンチマークビルド
├── run_all.sh                    # 全ベンチマーク実行
├── roofline_analysis.sh          # nsys/ncuルーフライン解析スクリプト
│
├── <id>_<kernel_name>/
│   ├── kernel_benchmark.f90      # 最適化版カーネル
│   ├── Makefile
│   ├── benchmark.conf
│   ├── data -> ../../Kernel_benchmark/<id>_.../data  # データはCPU版と共有
│   ├── OPTIMIZATION.md           # 最適化の記録（下記参照）
│   └── nsys_profile/             # nsysプロファイル結果
│
└── optimization_summary.csv      # 全カーネルの最適化結果一覧
```

---

## 作業手順

### Step 1: Kernel_benchmark_gpu_opt の作成

```bash
cp -r Kernel_benchmark_gpu Kernel_benchmark_gpu_opt
```

データディレクトリのシンボリックリンクが壊れていないか確認:

```bash
for d in Kernel_benchmark_gpu_opt/[0-9]*; do
  [ -L "$d/data" ] && [ ! -e "$d/data" ] && echo "BROKEN: $d/data"
done
```

### Step 2: nsys プロファイルの取得

各カーネルベンチマークに対して nsys プロファイルを取得する。

```bash
cd Kernel_benchmark_gpu_opt/<id>_<kernel_name>
nsys profile --stats=true -o nsys_profile/profile ./kernel_benchmark
```

#### nsys から得られる情報

| メトリクス | 意味 | 最適化への活用 |
|-----------|------|--------------|
| GPU kernel duration | カーネル実行時間 | ボトルネック特定 |
| Memory throughput | HBM帯域利用率 | メモリバウンドの判定 |
| SM occupancy | SM稼働率 | 並列度の充足判定 |
| Unified Memory page faults | UM転送イベント | ウォームアップの妥当性確認 |
| API call overhead | CUDA API呼び出しコスト | カーネル起動オーバーヘッド |

### Step 3: ルーフライン解析

既存の `roofline_analysis.sh`（ncu ベース）と nsys プロファイルを併用して各カーネルを分類:

| 分類 | Arithmetic Intensity | 最適化方針 |
|------|---------------------|-----------|
| メモリバウンド (AI < ridge point) | 低い | メモリアクセスパターン改善、帯域利用率向上 |
| 計算バウンド (AI >= ridge point) | 高い | 演算効率改善、不要計算の除去 |
| レイテンシバウンド | - | 並列度不足、occupancy 改善 |

### Step 4: 最適化の優先度決定

以下の基準で最適化対象カーネルを選定:

1. **実行時間の長いカーネル** — nsys のカーネル実行時間ランキング上位
2. **ピーク性能比が低いカーネル** — 理論ピークに対する達成率が低い
3. **最適化余地の大きいカーネル** — メモリアクセスが非効率、ループ構造が単純で改善しやすい

### Step 5: カーネル最適化の実施

`kernel_benchmark.f90` 内のカーネルサブルーチンのみを変更する。

#### 最適化テクニック一覧

**メモリバウンドカーネル向け:**

```fortran
! (A) collapse 節でループ並列度を向上
!$acc kernels
!$acc loop independent collapse(3)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      ...
    end do
  end do
end do
!$acc end kernels

! (B) 一時スカラで配列アクセス回数を削減
!$acc kernels
!$acc loop independent collapse(2)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      tmp_u = u(i,j,k)
      tmp_v = v(i,j,k)
      result(i,j,k) = tmp_u * tmp_v + tmp_u + tmp_v
    end do
  end do
end do
!$acc end kernels
```

**計算バウンドカーネル向け:**

```fortran
! (C) 除算を逆数乗算に置換
rdx = 1.0 / dx
!$acc kernels
!$acc loop independent collapse(3)
do k = 2, nk-2
  do j = 2, nj-2
    do i = 2, ni-2
      grad(i,j,k) = (phi(i+1,j,k) - phi(i,j,k)) * rdx
    end do
  end do
end do
!$acc end kernels
```

**Occupancy 改善:**

```fortran
! (D) tile 節でスレッドブロック形状を制御
!$acc kernels
!$acc loop independent tile(32, 4)
do j = 2, nj-2
  do i = 2, ni-2
    ...
  end do
end do
!$acc end kernels

! (E) gang/worker/vector の明示指定
!$acc parallel loop gang num_gangs(nk) vector_length(128)
do k = 2, nk-2
  !$acc loop vector
  do j = 2, nj-2
    do i = 2, ni-2
      ...
    end do
  end do
end do
```

### Step 6: 最適化の検証

最適化後、以下を確認:

1. **精度検証**: 既存の validation (reference 出力との比較) が PASSED であること
2. **性能比較**: 最適化前後の実行時間を比較

```bash
# 最適化前（gpu版）
cd Kernel_benchmark_gpu/<id>_<name> && ./kernel_benchmark

# 最適化後（gpu_opt版）
cd Kernel_benchmark_gpu_opt/<id>_<name> && ./kernel_benchmark
```

### Step 7: OPTIMIZATION.md の記録

各カーネルディレクトリに最適化内容を記録する:

```markdown
# Optimization Record: <kernel_name>

## Roofline Classification
- Arithmetic Intensity (DRAM): X.XX FLOP/Byte
- Classification: Memory-bound / Compute-bound
- Achieved BW: XXX GB/s (XX% of peak)
- Achieved GFLOP/s: XX.X (XX% of peak)

## Optimizations Applied
1. <最適化内容の簡潔な説明>
   - Before: <変更前のコードパターン>
   - After: <変更後のコードパターン>
   - Rationale: <最適化の根拠>

## Results
| Metric          | Before   | After    | Change  |
|-----------------|----------|----------|---------|
| Execution time  | X.XX ms  | X.XX ms  | -XX%    |
| Bandwidth       | XXX GB/s | XXX GB/s | +XX%    |
| Validation      | PASSED   | PASSED   | -       |
```

### Step 8: Src/ への統合

最適化が検証済みのカーネルを `Src/` の GPU ブランチ（`#if defined(USE_GPU)` 内）に反映する。

手順:
1. `Kernel_benchmark_gpu_opt/<id>/kernel_benchmark.f90` のカーネルサブルーチン部分を抽出
2. `Src/<source>.f90` の対応する `#if defined(USE_GPU)` ブロックを更新
3. `test_real/` で統合テスト実施（精度が維持されることを確認）

---

## optimization_summary.csv のフォーマット

```csv
kernel_id,kernel_name,classification,ai_dram,bw_before_gbs,bw_after_gbs,time_before_us,time_after_us,speedup,optimizations_applied,validation
067,copy3d,MEM,0.25,2100,2800,45.2,33.9,1.33x,"collapse(3)",PASSED
091,diver3d,MEM,0.41,1850,2600,62.1,44.3,1.40x,"collapse(3)+scalar_tmp",PASSED
...
```

---

## 注意事項

1. **1カーネルずつ最適化**: 複数カーネルを同時に変更しない。変更→検証→記録のサイクルを守る
2. **最適化前のコードを保持**: `Kernel_benchmark_gpu/` は変更しない（比較用ベースライン）
3. **入出力契約を守る**: サブルーチンの引数リスト（個数・型・形状・intent）を変更しない
4. **Unified Memory 前提を維持**: `!$acc data` 等のデータ配置ディレクティブを使わない
5. **精度劣化を許容しない**: validation が PASSED でなければその最適化はリバートする

---

## Related Documents

- [Phase 4: GPU Benchmark](04_gpu_benchmark.md) — 最適化前のベースライン
- [Phase 5: Simulation Integration](05_simulation_integration.md) — 統合テスト手順
- [Roofline Results](../../Kernel_benchmark_gpu/roofline_results/) — 既存のルーフライン解析結果

---

*Created: 2026-03-30*
