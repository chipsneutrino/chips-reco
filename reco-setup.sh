#! /bin/bash

if [ "x$CHIPSSIM" == "x" ]; then
    echo "You need to export CHIPSSIM to point to your chips-sim directory"
    return
fi

CURRENTDIR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export CHIPSRECO=$DIR

if [ -f "$DIR/evDisplay" ]
then
    echo "${C_GREEN}chips-reco built${C_RESET}"
else
    echo "${C_RED}building chips-reco${C_RESET}"
    NB_CORES=$(grep -c '^processor' /proc/cpuinfo)
    export MAKEFLAGS="-j$((NB_CORES+1)) -l${NB_CORES}"
    cd $DIR
    cmake .
    make
    cd $CURRENTDIR
    echo "${C_GREEN}chips-reco built${C_RESET}"
fi

export LD_LIBRARY_PATH=$DIR:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$DIR/include:$CPLUS_INCLUDE_PATH
export PATH=$DIR:$PATH
echo "${C_GREEN}chips-reco setup done${C_RESET}"