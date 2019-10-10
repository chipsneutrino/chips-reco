#! /bin/bash
if [ "x$WCSIMHOME" == "x" ]; then
    echo "You need to export WCSIMHOME to point to your WCSim directory"
    return
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export WCSIMANAHOME=${DIR}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/lib/
echo "WCSimAnalysis setup complete"
