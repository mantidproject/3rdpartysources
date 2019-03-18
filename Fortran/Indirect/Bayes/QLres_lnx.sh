#!/bin/sh

f2py -c --fcompiler=gnu95 --compiler=unix \
 -m QLres_lnx64 \
  QLres_main.f90 only: qlres : QLres_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90
