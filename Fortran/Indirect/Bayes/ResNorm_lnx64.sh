#!/bin/sh

f2py3 -c --fcompiler=gnu95 --compiler=unix \
 -m ResNorm_lnx64 \
  ResNorm_main.f90 only: resnorm : ResNorm_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90
