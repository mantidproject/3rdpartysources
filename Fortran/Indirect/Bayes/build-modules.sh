#!/bin/bash
# Build all fortran Python modules

LOGSDIR="build_logs"
mkdir -p $LOGSDIR

function exit_if_failed() {
  modulename=$1
  result=$2
  if [ $2 -ne 0 ]; then
    echo "Building $modulename failed. See build log in $LOGSDIR"
    exit 1
  fi
}

function build_module() {
  modulename=$1
  mainsrc=$2
  only=$3
  othersrcs=${@:5}
  echo "Building $modulename"
  f2py3 -c --fcompiler=gnu95 --compiler=unix \
    -m $modulename \
    $mainsrc only: $only : $othersrcs > $LOGSDIR/$modulename.log 2>&1
  exit_if_failed $modulename $?
}


# QLdata
build_module QLdata QLdata_main.f90 qldata QLdata_subs.f90 Bayes.f90 Four.f90 Util.f90
# QLres
build_module QLres QLres_main.f90 qlres QLres_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90
# QLse
build_module QLse QLse_main.f90 qlstexp QLse_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90 Simopt.f90
# Quest
build_module Quest Quest_main.f90 quest Quest_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90 Simopt.f90
# ResNorm
build_module ResNorm ResNorm_main.f90 resnorm ResNorm_subs.f90 BlrRes.f90 Bayes.f90 Four.f90 Util.f90

echo "All modules built successfully"
