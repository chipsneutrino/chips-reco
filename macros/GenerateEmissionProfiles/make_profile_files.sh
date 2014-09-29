#!/bin/bash
##########################
#   the script to run   # .mac name # .root name # particle # energy # nevents # split files after nevts # WCSim directory
python generate_macfile.py -o `pwd`/out/ep_1000 -r `pwd`/out/ep_muon_1000 -p mu- -e 1000 -n 100 -s 100 -w ${WCSIMHOME}
python generate_macfile.py -o `pwd`/out/ep_1500 -r `pwd`/out/ep_muon_1500 -p mu- -e 1500 -n 100 -s 100 -w ${WCSIMHOME}
python generate_macfile.py -o `pwd`/out/ep_2000 -r `pwd`/out/ep_muon_2000 -p mu- -e 2000 -n 100 -s 100 -w ${WCSIMHOME}
python generate_macfile.py -o `pwd`/out/ep_2500 -r `pwd`/out/ep_muon_2500 -p mu- -e 2500 -n 100 -s 100 -w ${WCSIMHOME}
python generate_macfile.py -o `pwd`/out/ep_3000 -r `pwd`/out/ep_muon_3000 -p mu- -e 3000 -n 100 -s 100 -w ${WCSIMHOME}
