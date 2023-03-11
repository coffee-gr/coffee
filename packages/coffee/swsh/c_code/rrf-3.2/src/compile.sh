#!/bin/bash

FC=$1
FFLAGS=""
shift
while [[ $# -gt 0 ]]; do
  set FFLAGS = "$FFLAGS $1"
  shift
done

if [[ $FC = "f95" || $FC = "nagfor" ]]; then
  DEFS="-DNAGF95"
elif [[ $FC = g95 ]]; then
  DEFS="-DG95"
fi

$FC $FFLAGS $DEFS -c rrfdata.F90
$FC $FFLAGS -c rrflist.f90
$FC $FFLAGS -c arithmetic.f90
$FC $FFLAGS -c factorials.f90
$FC $FFLAGS -c start.f90
$FC $FFLAGS -c conversions.f90
$FC $FFLAGS -c wigner.F90
$FC $FFLAGS $DEFS -c input.F90
