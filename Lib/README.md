# Kernel I/O Library

カーネルベンチマーク用の統一されたI/Oライブラリです。

## 概要

このライブラリは以下の機能を提供します：

1. **データ書き出し (Dump)**: シミュレーション実行中にカーネルの入出力データをバイナリファイルに保存
2. **データ読み込み (Load)**: ベンチマークプログラムでダンプされたデータを読み込み
3. **1GB バッファ制限チェック**: Fortran I/O の内部バッファ制限（約1GB）を超える配列に対するエラー検出

## バッファ制限と自動分割

### 背景

Fortran の一部のランタイム実装では、配列 I/O に内部バッファ制限があります。
この制限はコンパイラや環境によって異なります：

| コンパイラ | 状況 |
|------------|------|
| gfortran | 環境変数 `GFORTRAN_UNBUFFERED_ALL` 等で制御可能 |
| Intel Fortran | `FORT_BUFFERED`, `FORT_BLOCKSIZE` 等で制御可能 |
| nvfortran | 独自の制限がある場合あり |

### デフォルト制限値

デフォルトでは保守的に **1GB** を制限値として使用します：

- 4バイトデータ (real, integer): 最大 268,435,456 要素
- 8バイトデータ (real(8), integer(8)): 最大 134,217,728 要素

### 制限値のカスタマイズ

環境や用途に応じて制限値を変更できます：

**方法1: 環境変数**
```bash
export KIO_BUFFER_LIMIT_MB=512   # 512MBに設定
./your_program
```

**方法2: プログラム内で設定**
```fortran
use m_kernel_io
call kio_set_buffer_limit(512)   ! 512MBに設定
```

**方法3: 環境変数を明示的に読み込み**
```fortran
use m_kernel_io
logical :: found
found = kio_init_from_env()      ! KIO_BUFFER_LIMIT_MB を読み込み
```

### 自動分割機能

配列サイズが制限を超える場合、自動的にブロック分割してI/Oを行います：

| 配列次元 | 分割方式 |
|----------|----------|
| 1D | チャンク単位で分割 |
| 2D | 1Dスライス (j方向ループ) |
| 3D | 2Dスライス (k方向ループ) |
| 4D | 常に2Dスライス (k,l方向ループ) |

分割が行われた場合、ログに表示されます：
```
[KIO] Buffer limit set to 512 MB
[KIO] 3D: large_array.bin (    500000000) [2D slices]
```

**注意**: 2Dスライス自体が制限を超える場合（i×j が制限要素数以上）はエラーとなります。

## 使用方法

### ビルド

```bash
cd Lib
make          # ライブラリをビルド
make test     # テストを実行
make clean    # クリーンアップ
```

### データ書き出し（Dump）

```fortran
use m_kernel_io

! 初期化
call kio_dump_init('kernel_name')

! スカラー値の保存
call kio_dump_scalar_i('ni', ni)
call kio_dump_scalar_r('dx', dx)
call kio_dump_scalar_d('pi', 3.14159265358979d0)
call kio_dump_scalar_c('name', 'value')

! 配列の保存
call kio_dump_array_1d('arr1d.bin', arr, 1, n)
call kio_dump_array_2d('arr2d.bin', arr, 0, ni+1, 0, nj+1)
call kio_dump_array_3d('arr3d.bin', arr, 0, ni+1, 0, nj+1, 1, nk)
call kio_dump_array_4d('arr4d.bin', arr, 0, ni+1, 0, nj+1, 1, nk, 1, nq)

! 整数配列
call kio_dump_array_1d_int('iarr.bin', iarr, 1, n)

! 終了
call kio_dump_finalize()
```

### データ読み込み（Load）

```fortran
use m_kernel_io

! パラメータファイルの読み込み
call kio_load_params('data/params.txt')

! パラメータの取得
call kio_get_param_i('ni', ni)
call kio_get_param_r('dx', dx)
call kio_get_param_d('pi', pi)
call kio_get_param_c('name', name)

! 配列の読み込み
call kio_load_array_1d('data/arr1d.bin', arr, 1, n)
call kio_load_array_2d('data/arr2d.bin', arr, 0, ni+1, 0, nj+1)
call kio_load_array_3d('data/arr3d.bin', arr, 0, ni+1, 0, nj+1, 1, nk)
call kio_load_array_4d('data/arr4d.bin', arr, 0, ni+1, 0, nj+1, 1, nk, 1, nq)
```

### 後方互換性

既存の `m_dump_kernel` モジュールを使用しているコードは、そのまま動作します：

```fortran
use m_dump_kernel  ! 従来通り動作

call dump_init('kernel_name')
call dump_scalar_i('ni', ni)
call dump_array_3d('u.bin', u, 0, ni+1, 0, nj+1, 1, nk)
call dump_finalize()
```

### サイズチェック関数

プログラム内で明示的にサイズをチェックすることもできます：

```fortran
use m_kernel_io

integer(8) :: nelements

nelements = int(ni+2, 8) * int(nj+2, 8) * int(nk, 8)

if (.not. kio_check_size(nelements)) then
  write(*,*) 'Warning: Array exceeds 1GB buffer limit!'
  ! エラー処理
end if
```

## API リファレンス

### 定数

| 名前 | 値 | 説明 |
|------|-----|------|
| `KIO_BUFFER_LIMIT_BYTES` | 1,073,741,824 | 1GB (バイト単位) |
| `KIO_MAX_ELEMENTS_4BYTE` | 268,435,456 | 4バイト要素の最大数 |
| `KIO_MAX_ELEMENTS_8BYTE` | 134,217,728 | 8バイト要素の最大数 |

### Dump ルーチン

| ルーチン | 説明 |
|----------|------|
| `kio_dump_init(kernel_name)` | ダンプ初期化 |
| `kio_dump_finalize()` | ダンプ終了 |
| `kio_dump_scalar_i(name, val)` | 整数スカラー |
| `kio_dump_scalar_i8(name, val)` | 8バイト整数スカラー |
| `kio_dump_scalar_r(name, val)` | 単精度実数スカラー |
| `kio_dump_scalar_d(name, val)` | 倍精度実数スカラー |
| `kio_dump_scalar_c(name, val)` | 文字列スカラー |
| `kio_dump_array_1d(file, arr, i1, i2)` | 1D実数配列 |
| `kio_dump_array_2d(file, arr, i1, i2, j1, j2)` | 2D実数配列 |
| `kio_dump_array_3d(file, arr, i1, i2, j1, j2, k1, k2)` | 3D実数配列 |
| `kio_dump_array_4d(file, arr, ...)` | 4D実数配列（2Dスライス） |
| `kio_dump_array_Nd_int(...)` | 整数配列 (N=1,2,3) |

### Load ルーチン

| ルーチン | 説明 |
|----------|------|
| `kio_load_params(filename)` | パラメータファイル読み込み |
| `kio_params_loaded()` | 読み込み済みか確認 |
| `kio_get_param_i(name, val [, default])` | 整数パラメータ取得 |
| `kio_get_param_r(name, val [, default])` | 単精度実数パラメータ取得 |
| `kio_get_param_d(name, val [, default])` | 倍精度実数パラメータ取得 |
| `kio_get_param_c(name, val [, default])` | 文字列パラメータ取得 |
| `kio_load_array_Nd(file, arr, bounds...)` | 配列読み込み (N=1,2,3,4) |

### ユーティリティ

| ルーチン | 説明 |
|----------|------|
| `kio_check_size(nelements)` | 要素数が制限内かチェック（4バイト用） |
| `kio_check_size_bytes(nbytes)` | バイト数が制限内かチェック |
| `kio_error(msg)` | エラーメッセージを表示して停止 |

### バッファ制限設定

| ルーチン | 説明 |
|----------|------|
| `kio_set_buffer_limit(limit_mb)` | バッファ制限をMB単位で設定 |
| `kio_get_buffer_limit()` | 現在のバッファ制限（バイト）を取得 |
| `kio_init_from_env()` | 環境変数 `KIO_BUFFER_LIMIT_MB` から初期化 |

## ファイル構成

```
Lib/
├── m_kernel_io.f90      # メインライブラリモジュール
├── m_dump_kernel.f90    # 後方互換性ラッパー
├── test_kernel_io.f90   # テストプログラム
├── Makefile             # ビルド用 Makefile
└── README.md            # このファイル
```

## エラーメッセージ例

配列サイズが1GB制限を超えた場合：

```
***** KERNEL I/O ERROR *****
Array size exceeds 1GB I/O buffer limit!
  Array: large_array.bin
  Dimensions: (0:1000, 0:1000, 1:500)
  Elements: 501501500 (max: 268435456)
  Bytes: 2006006000 (max: 1073741824)
****************************
```
