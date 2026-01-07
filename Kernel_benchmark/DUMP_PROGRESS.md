# Kernel Data Dump Progress

## Last Updated: 2026-01-03

## Summary
Srcフォルダのソースファイルに、最終呼び出し時のデータダンプコードを追加完了

## dump_kernel_data.f90
- 場所: Src/dump_kernel_data.f90
- dump_scalar_c (文字列用) を追加済み

## 完了したダンプコード追加 (全20個完了)
| ファイル | 実行回数 | ダンプ対象呼び出し | ステータス |
|----------|----------|-------------------|-----------|
| advbspi.f90 | 28800 | 28800 | 完了 |
| diver2d.f90 | 14400 | 14400 | 完了 |
| advs.f90 | 3960 | 3960 | 完了 |
| buoywsi.f90 | 14400 | 14400 | 完了 |
| diver3d.f90 | 14400 | 14400 | 完了 |
| pgrad.f90 | 14400 | 14400 | 完了 |
| stepuv.f90 | 14400 | 14400 | 完了 |
| stepwi.f90 | 14400 | 14400 | 完了 |
| steppi.f90 | 14400 | 14400 | 完了 |
| phy2cnt.f90 | 15121 | 15121 | 完了 |
| vbcu.f90 | 14401 | 14401 | 完了 |
| vbcv.f90 | 14401 | 14401 | 完了 |
| vbcw.f90 | 14401 | 14401 | 完了 |
| vbcwc.f90 | 15121 | 15121 | 完了 |
| steps.f90 | 360 | 360 | 完了 |
| disptke.f90 | 360 | 360 | 完了 |
| copy3d.f90 | 1446 | 1446 | 完了 |
| bcycle.f90 | 72374 | 72374 | 完了 |
| bc4news.f90 | 72374 | 72374 | 完了 |
| diagni.f90 | 1 | 1 | 完了 |

## ダンプコード追加パターン

```fortran
! モジュール参照に追加
use m_dump_kernel

! 変数宣言部に追加
integer, save :: dump_call_count_XXXX = 0
integer, parameter :: DUMP_TARGET_XXXX = <count>
logical, save :: dump_done_XXXX = .false.

! profile_start直前に追加
dump_call_count_XXXX = dump_call_count_XXXX + 1
if (dump_call_count_XXXX == DUMP_TARGET_XXXX .and. .not. dump_done_XXXX) then
  call dump_init('XXXX')
  call dump_scalar_i('ni', ni)
  ! ... 入力データのダンプ
end if

! profile_stop直後に追加
if (dump_call_count_XXXX == DUMP_TARGET_XXXX .and. .not. dump_done_XXXX) then
  call dump_array_3d('output_ref.bin', output, ...)
  call dump_finalize()
  dump_done_XXXX = .true.
end if
```

## ダンプデータ生成手順
1. ダンプコード追加済みのSrcでシミュレーションをビルド
2. test_realで実行
3. ./kernel_dump/<kernel_name>/ にデータが生成される
4. 各ベンチマークのdata/にコピー

## 次のステップ
1. 親ディレクトリでmakeを実行してビルド
2. test_realフォルダでシミュレーションを実行
3. kernel_dumpディレクトリに生成されたデータを確認
4. 各カーネルベンチマークのdata/ディレクトリにコピー
