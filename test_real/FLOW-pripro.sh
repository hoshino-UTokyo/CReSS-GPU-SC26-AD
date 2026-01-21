#!/bin/sh
# this script can be only used FX10.
#PJM --rsc-list "rscunit=fx"
#PJM --rsc-list "rscgrp=fx-small"
#PJM --rsc-list "node=1,elapse=10:00:00,node-mem=25600"
#PJM --name "pripro.sh"
#PJM -j
#PJM -S

#./check.exe < user.conf > check.txt

rm -f result/*gpv*
./gridata.exe -Wl,-Lu -Wl,-T < user.conf > log.gridata.txt
