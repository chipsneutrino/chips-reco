#! /bin/bash

if [ "x$WCSIMHOME" == "x" ]; then
    echo "You need to export WCSIMHOME to point to your WCSim directory"
    return
fi


#WCSimAnalysis Home Directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export WCSIMANAHOME=${DIR}

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/lib/

#FIXME: necessary? just run $WCS/setupWCSim.sh?
source /unix/lartpc/software/root/setup.sh
echo "Root setup complete"
source /unix/lartpc/software/geant4/setup.sh
echo "geant4 setup complete"
source /unix/lartpc/software/genie/setup.sh
echo "genie setup complete"
