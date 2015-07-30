#! /bin/bash

export WCSIMHOME=/unix/fnu/ajperch/software/WCSim_github

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
