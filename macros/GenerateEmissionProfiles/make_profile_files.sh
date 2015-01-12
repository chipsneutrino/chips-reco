#!/bin/bash
##########################
#   the script to run   # .mac name # .root name # particle # energy # nevents # split files after nevts # WCSim directory
python generate_macfile.py -o `pwd`/out/ep_1250_fast -r ${FNU}/ep_jan5/ep_muon_1250_fast -p mu- -e 1250 -n 5 -s 5 -w ${WCSIMHOME}
python generate_macfile.py -o `pwd`/out/ep_1500_fast -r ${FNU}/ep_jan5/ep_muon_1500_fast -p mu- -e 1500 -n 5 -s 5 -w ${WCSIMHOME}
python generate_macfile.py -o `pwd`/out/ep_1750_fast -r ${FNU}/ep_jan5/ep_muon_1750_fast -p mu- -e 1750 -n 5 -s 5 -w ${WCSIMHOME}
