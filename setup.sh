#! /bin/bash

if [ "x$WCSIMHOME" == "x" ]; then
    echo "You need to export WCSIMHOME to point to your WCSim directory"
    return
fi

# Set the chips-reco directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export WCSIMANAHOME=$DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/lib/
export PATH=$PATH:$DIR/