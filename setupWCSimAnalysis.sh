#! /bin/bash

#WCSimAnalysis Home Directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export WCSIMANAHOME=${DIR}

#FIXME: necessary? just run $WCS/setupWCSim.sh?
source /unix/lartpc/software/root/setup.sh
echo "Root setup complete"
source /unix/lartpc/software/geant4/setup.sh
echo "geant4 setup complete"
source /unix/lartpc/software/genie/setup.sh
echo "genie setup complete"
