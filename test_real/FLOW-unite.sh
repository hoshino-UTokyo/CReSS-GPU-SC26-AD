#!/bin/sh
# this script can be only used FX10.
#PJM --rsc-list "rscunit=fx"
#PJM --rsc-list "rscgrp=fx-small"
#PJM --rsc-list "node=1,elapse=10:00:00,node-mem=25600"
#PJM --name "log-unite.sh"
#PJM -j
#PJM -S

#./check.exe < user.conf > check.txt

./unite.exe -Wl,-Lu -Wl,-T < user.conf > log.unite.txt

#./unite.exe -Wl,-Lu -Wl,-T < user.conf_geo > log.unite-geo.txt
#./unite.exe -Wl,-Lu -Wl,-T < user.conf_mon > log.unite-mon2.txt
#./unite.exe -Wl,-Lu -Wl,-T < user.conf_dmp > log.unite-dmp.txt
