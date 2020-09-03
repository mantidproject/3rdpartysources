#!/bin/sh

f2py3 -c --fcompiler=gnu95 --compiler=unix \
 -m QLdata_lnx64 \
  QLdata_main.f90 only: qldata : QLdata_subs.f90 Bayes.f90 Four.f90 Util.f90
