#!/bin/sh
cd $DIR

for (( BIN=0; BIN<${BINSPERJOB}; BIN++ ))
do
  WHICHBIN=$((MUBIN+BIN))
  ARG="(${WHICHBIN}, ${NUMTHROWS}, ${NUMBINS}, ${MIN}, ${MAX})"
  root -b -q "RunDigitizerPDFMaker.C$ARG"
done
