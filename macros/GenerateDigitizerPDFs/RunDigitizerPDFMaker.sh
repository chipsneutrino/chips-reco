#!/bin/sh
cd $DIR

source /unix/chips/jtingey/CHIPS/code/WCSim/setupScripts/setupWCSimUCL.sh
source /unix/chips/jtingey/CHIPS/code/WCSimAnalysis/setupWCSimAnalysis.sh
for (( BIN=0; BIN<${BINSPERJOB}; BIN++ ))
do
  WHICHBIN=$((MUBIN+BIN))
  ARG="(${WHICHBIN}, ${TYPE}, ${NUMTHROWS}, ${NUMBINS}, ${MIN}, ${MAX})"
  root -l -b -q "RunDigitizerPDFMaker.C$ARG"
done
