# Kernel Dump Data Summary

This document provides a quick reference for all kernel dump configurations.

## Overview

- **Total Kernels with Dump Code**: 157 (out of 387 kernel directories)
- **Dump Directory Base**: `test_real/kernel_dump/`
- **Dump Configuration**: Each kernel dumps at its final call (from omp_profile.txt)

## Quick Reference Table

| Kernel | Source File | Dump Target | Dump Dir |
|--------|------------|-------------|----------|
| adjstni | adjstni.f90 | 1080 | `kernel_dump/adjstni/` |
| adjstq | adjstq.f90 | 3 | `kernel_dump/adjstq/` |
| adjstuv | adjstuv.f90 | 360 | `kernel_dump/adjstuv/` |
| advbspi | advbspi.f90 | 28800 | `kernel_dump/advbspi/` |
| advbspt | advbspt.f90 | 360 | `kernel_dump/advbspt/` |
| advp | advp.f90 | 360 | `kernel_dump/advp/` |
| advs | advs.f90 | 3960 | `kernel_dump/advs/` |
| advuvw | advuvw.f90 | 360 | `kernel_dump/advuvw/` |
| aggregat | aggregat.f90 | 45720 | `kernel_dump/aggregat/` |
| allocbuf | allocbuf.f90 | 1 | `kernel_dump/allocbuf/` |
| allociot | allociot.f90 | 1 | `kernel_dump/allociot/` |
| allocslv | allocslv.f90 | 1 | `kernel_dump/allocslv/` |
| baserho | baserho.f90 | 1 | `kernel_dump/baserho/` |
| bbcw | bbcw.f90 | 14400 | `kernel_dump/bbcw/` |
| bc4news | bc4news.f90 | 72374 | `kernel_dump/bc4news/` |
| bc8u | bc8u.f90 | 720 | `kernel_dump/bc8u/` |
| bc8v | bc8v.f90 | 720 | `kernel_dump/bc8v/` |
| bc8w | bc8w.f90 | 2 | `kernel_dump/bc8w/` |
| bcbase | bcbase.f90 | 1 | `kernel_dump/bcbase/` |
| bcten | bcten.f90 | 1440 | `kernel_dump/bcten/` |
| bcycle | bcycle.f90 | 72374 | `kernel_dump/bcycle/` |
| bcyclex | bcyclex.f90 | 4 | `kernel_dump/bcyclex/` |
| bcycley | bcycley.f90 | 4 | `kernel_dump/bcycley/` |
| bruntv | bruntv.f90 | 360 | `kernel_dump/bruntv/` |
| bulksfc | bulksfc.f90 | 385 | `kernel_dump/bulksfc/` |
| buoytke | buoytke.f90 | 360 | `kernel_dump/buoytke/` |
| buoywb | buoywb.f90 | 360 | `kernel_dump/buoywb/` |
| buoywsi | buoywsi.f90 | 14400 | `kernel_dump/buoywsi/` |
| chkitr | chkitr.f90 | 16 | `kernel_dump/chkitr/` |
| chkmoist | chkmoist.f90 | 1 | `kernel_dump/chkmoist/` |
| chkrain | chkrain.f90 | 361 | `kernel_dump/chkrain/` |
| chksat | chksat.f90 | 3 | `kernel_dump/chksat/` |
| cloudcov | cloudcov.f90 | 361 | `kernel_dump/cloudcov/` |
| collect | collect.f90 | 45720 | `kernel_dump/collect/` |
| convers | convers.f90 | 45720 | `kernel_dump/convers/` |
| copy1d | copy1d.f90 | 1 | `kernel_dump/copy1d/` |
| copy3d | copy3d.f90 | 1446 | `kernel_dump/copy3d/` |
| copy4d | copy4d.f90 | 3 | `kernel_dump/copy4d/` |
| coriuv | coriuv.f90 | 360 | `kernel_dump/coriuv/` |
| curveuvw | curveuvw.f90 | 360 | `kernel_dump/curveuvw/` |
| defomssq | defomssq.f90 | 360 | `kernel_dump/defomssq/` |
| depsit | depsit.f90 | 45720 | `kernel_dump/depsit/` |
| diagni | diagni.f90 | 1 | `kernel_dump/diagni/` |
| diagnw | diagnw.f90 | 360 | `kernel_dump/diagnw/` |
| disptke | disptke.f90 | 360 | `kernel_dump/disptke/` |
| distrpg | distrpg.f90 | 45720 | `kernel_dump/distrpg/` |
| diver2d | diver2d.f90 | 14400 | `kernel_dump/diver2d/` |
| diver3d | diver3d.f90 | 14400 | `kernel_dump/diver3d/` |
| diverpih | diverpih.f90 | 14400 | `kernel_dump/diverpih/` |
| diverpiv | diverpiv.f90 | 28800 | `kernel_dump/diverpiv/` |
| eddydif | eddydif.f90 | 360 | `kernel_dump/eddydif/` |
| eddyvis | eddyvis.f90 | 360 | `kernel_dump/eddyvis/` |
| eddyvisj | eddyvisj.f90 | 360 | `kernel_dump/eddyvisj/` |
| exbcpt | exbcpt.f90 | 360 | `kernel_dump/exbcpt/` |
| exbcq | exbcq.f90 | 360 | `kernel_dump/exbcq/` |
| exbcss | exbcss.f90 | 14400 | `kernel_dump/exbcss/` |
| exbcu | exbcu.f90 | 14400 | `kernel_dump/exbcu/` |
| exbcv | exbcv.f90 | 14400 | `kernel_dump/exbcv/` |
| fallblk | fallblk.f90 | 360 | `kernel_dump/fallblk/` |
| forcept | forcept.f90 | 360 | `kernel_dump/forcept/` |
| forcesfc | forcesfc.f90 | 361 | `kernel_dump/forcesfc/` |
| freezing | freezing.f90 | 45720 | `kernel_dump/freezing/` |
| gaussel | gaussel.f90 | 14418 | `kernel_dump/gaussel/` |
| getarea | getarea.f90 | 1 | `kernel_dump/getarea/` |
| getexner | getexner.f90 | 1080 | `kernel_dump/getexner/` |
| getrich | getrich.f90 | 377 | `kernel_dump/getrich/` |
| gettrn | gettrn.f90 | 1 | `kernel_dump/gettrn/` |
| getvdens | getvdens.f90 | 360 | `kernel_dump/getvdens/` |
| getxy | getxy.f90 | 2 | `kernel_dump/getxy/` |
| getz | getz.f90 | 1 | `kernel_dump/getz/` |
| getzlow | getzlow.f90 | 361 | `kernel_dump/getzlow/` |
| heatsfc | heatsfc.f90 | 361 | `kernel_dump/heatsfc/` |
| inidef | inidef.f90 | 1 | `kernel_dump/inidef/` |
| inimod | inimod.f90 | 1 | `kernel_dump/inimod/` |
| inisfc | inisfc.f90 | 1 | `kernel_dump/inisfc/` |
| initund | initund.f90 | 1 | `kernel_dump/initund/` |
| jacobian | jacobian.f90 | 1 | `kernel_dump/jacobian/` |
| kh8uv | kh8uv.f90 | 360 | `kernel_dump/kh8uv/` |
| lbcs | lbcs.f90 | 3240 | `kernel_dump/lbcs/` |
| lbcw | lbcw.f90 | 14400 | `kernel_dump/lbcw/` |
| mapfct | mapfct.f90 | 1 | `kernel_dump/mapfct/` |
| melting | melting.f90 | 45720 | `kernel_dump/melting/` |
| more0q | more0q.f90 | 45720 | `kernel_dump/more0q/` |
| newblk | newblk.f90 | 45720 | `kernel_dump/newblk/` |
| nuc1stc | nuc1stc.f90 | 45720 | `kernel_dump/nuc1stc/` |
| nuc1stv | nuc1stv.f90 | 45720 | `kernel_dump/nuc1stv/` |
| nuc2nd | nuc2nd.f90 | 45720 | `kernel_dump/nuc2nd/` |
| opendmp | opendmp.f90 | 4 | `kernel_dump/opendmp/` |
| outdmp | outdmp.f90 | 4 | `kernel_dump/outdmp/` |
| outmxn | outmxn.f90 | 5415 | `kernel_dump/outmxn/` |
| outpbl | outpbl.f90 | 4 | `kernel_dump/outpbl/` |
| pgrad | pgrad.f90 | 14400 | `kernel_dump/pgrad/` |
| pgradiv | pgradiv.f90 | 14400 | `kernel_dump/pgradiv/` |
| phvbcs | phvbcs.f90 | 1080 | `kernel_dump/phvbcs/` |
| phvbcuvw | phvbcuvw.f90 | 360 | `kernel_dump/phvbcuvw/` |
| phvs | phvs.f90 | 4320 | `kernel_dump/phvs/` |
| phvuvw | phvuvw.f90 | 360 | `kernel_dump/phvuvw/` |
| phy2cnt | phy2cnt.f90 | 15121 | `kernel_dump/phy2cnt/` |
| phycood | phycood.f90 | 1 | `kernel_dump/phycood/` |
| prodctwg | prodctwg.f90 | 45720 | `kernel_dump/prodctwg/` |
| radiat | radiat.f90 | 361 | `kernel_dump/radiat/` |
| rbcq | rbcq.f90 | 1800 | `kernel_dump/rbcq/` |
| rbcs | rbcs.f90 | 1080 | `kernel_dump/rbcs/` |
| rbcs0 | rbcs0.f90 | 360 | `kernel_dump/rbcs0/` |
| rbcw | rbcw.f90 | 14400 | `kernel_dump/rbcw/` |
| rdgrp | rdgrp.f90 | 1 | `kernel_dump/rdgrp/` |
| roughitr | roughitr.f90 | 16 | `kernel_dump/roughitr/` |
| roughnxt | roughnxt.f90 | 361 | `kernel_dump/roughnxt/` |
| rstuvwc | rstuvwc.f90 | 360 | `kernel_dump/rstuvwc/` |
| setbase | setbase.f90 | 1 | `kernel_dump/setbase/` |
| setblk | setblk.f90 | 45720 | `kernel_dump/setblk/` |
| setcst1d | setcst1d.f90 | 1 | `kernel_dump/setcst1d/` |
| setcst2d | setcst2d.f90 | 41 | `kernel_dump/setcst2d/` |
| setcst3d | setcst3d.f90 | 108 | `kernel_dump/setcst3d/` |
| setcst4d | setcst4d.f90 | 25 | `kernel_dump/setcst4d/` |
| setgpv | setgpv.f90 | 2 | `kernel_dump/setgpv/` |
| setname | setname.f90 | 1 | `kernel_dump/setname/` |
| setsfc | setsfc.f90 | 361 | `kernel_dump/setsfc/` |
| sfcflx | sfcflx.f90 | 361 | `kernel_dump/sfcflx/` |
| sheartke | sheartke.f90 | 360 | `kernel_dump/sheartke/` |
| shedding | shedding.f90 | 45720 | `kernel_dump/shedding/` |
| siadjst | siadjst.f90 | 720 | `kernel_dump/siadjst/` |
| smoo4qv | smoo4qv.f90 | 360 | `kernel_dump/smoo4qv/` |
| smoo4s | smoo4s.f90 | 3600 | `kernel_dump/smoo4s/` |
| smoo4uvw | smoo4uvw.f90 | 360 | `kernel_dump/smoo4uvw/` |
| sndwave | sndwave.f90 | 1 | `kernel_dump/sndwave/` |
| steppi | steppi.f90 | 14400 | `kernel_dump/steppi/` |
| steps | steps.f90 | 360 | `kernel_dump/steps/` |
| steptund | steptund.f90 | 18 | `kernel_dump/steptund/` |
| stepuv | stepuv.f90 | 14400 | `kernel_dump/stepuv/` |
| stepwi | stepwi.f90 | 14400 | `kernel_dump/stepwi/` |
| strsten | strsten.f90 | 360 | `kernel_dump/strsten/` |
| swadjst | swadjst.f90 | 720 | `kernel_dump/swadjst/` |
| swp2nxt | swp2nxt.f90 | 360 | `kernel_dump/swp2nxt/` |
| termblk | termblk.f90 | 360 | `kernel_dump/termblk/` |
| timeflt | timeflt.f90 | 359 | `kernel_dump/timeflt/` |
| totalqwi | totalqwi.f90 | 720 | `kernel_dump/totalqwi/` |
| totals | totals.f90 | 8 | `kernel_dump/totals/` |
| trilat | trilat.f90 | 1 | `kernel_dump/trilat/` |
| turbs | turbs.f90 | 3600 | `kernel_dump/turbs/` |
| turbtke | turbtke.f90 | 360 | `kernel_dump/turbtke/` |
| turbuvw | turbuvw.f90 | 360 | `kernel_dump/turbuvw/` |
| upwnp | upwnp.f90 | 1080 | `kernel_dump/upwnp/` |
| upwqp | upwqp.f90 | 1800 | `kernel_dump/upwqp/` |
| var8uvw | var8uvw.f90 | 2 | `kernel_dump/var8uvw/` |
| vbcp | vbcp.f90 | 14401 | `kernel_dump/vbcp/` |
| vbcs | vbcs.f90 | 3962 | `kernel_dump/vbcs/` |
| vbcu | vbcu.f90 | 14401 | `kernel_dump/vbcu/` |
| vbcv | vbcv.f90 | 14401 | `kernel_dump/vbcv/` |
| vbcw | vbcw.f90 | 14401 | `kernel_dump/vbcw/` |
| vbcwc | vbcwc.f90 | 15121 | `kernel_dump/vbcwc/` |
| vspdmp | vspdmp.f90 | 1 | `kernel_dump/vspdmp/` |
| vspqv | vspqv.f90 | 360 | `kernel_dump/vspqv/` |
| vsps | vsps.f90 | 2160 | `kernel_dump/vsps/` |
| vsps0 | vsps0.f90 | 1440 | `kernel_dump/vsps0/` |
| vspuvw | vspuvw.f90 | 360 | `kernel_dump/vspuvw/` |
| xy2ll | xy2ll.f90 | 1 | `kernel_dump/xy2ll/` |

## How Dump Works

1. **Dump Trigger**: Each kernel dumps data on its final call (DUMP_TARGET)
2. **Dump Location**: Data is saved to `test_real/kernel_dump/<kernel_name>/`
3. **File Types**:
   - `params.txt`: Scalar parameters (grid sizes, options, constants)
   - `*.bin`: Binary array files (Fortran stream format)
   - `*_ref.bin`: Reference output for validation

## Usage

### Running Simulation with Dump
```bash
cd test_real
./run_solver.sh  # Dumps will be created automatically
```

### Checking Dump Status
```bash
# List all dumped kernels
ls test_real/kernel_dump/

# Check specific kernel dump
ls -la test_real/kernel_dump/advbspi/
```

### Linking Data to Benchmark
```bash
# Create symlink in benchmark directory
cd Kernel_benchmark/012_advbspi_subroutine/
ln -s ../../test_real/kernel_dump/advbspi data
```

## Notes

- Dump code is added to **157 source files** covering all active kernels
- Each kernel's README.md contains detailed dump requirements
- Some kernels in omp_profile with Count=0 are not executed in test_real and won't produce dumps
