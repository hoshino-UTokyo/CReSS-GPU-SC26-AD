#!/bin/bash
################################################################################
#                                                                              #
#     ISPACK download script.                                                  #
#                                                                              #
#     Author      : Satoki Tsujino                                             #
#     Date        : 2017/06/14                                                 #
#     Modification:                                                            #
#                                                                              #
################################################################################

dir='RAD1.4.0'

if [ $1 = 'download' ]; then 
  wget http://www.gfd-dennou.org/library/ispack/ispack-1.0.4.tar.gz
  tar zxvf ispack-1.0.4.tar.gz

  for i in appack bspack c2pack fepack fhpack flpack ftpack p2pack 
  do
    cp "$i"/src/*.{f,c,F} SrcDevelop/"$dir"/10_fft_ispack/
  done
  rm -rf ispack-1.0.4 ispack-1.0.4.tar.gz
fi

if [ $1 = 'clean' ]; then
  rm -rf SrcDevelop/"$dir"/10_fft_ispack/*
  echo "rm -rf SrcDevelop/$dir/10_fft_ispack/*"
fi
