#!/usr/bin/python
import sys
import getopt
import os
from math import ceil


def main():
  
  ###############################
  # o == option
  # a == argument passed to the o
  ###############################
  # Cache an error with try..except
  # Note: options is the string of option letters that the script wants to recognize, with
  # options that require an argument followed by a colon (':') i.e. -i fileName# Short options are one letter, longer equivalents in the square brackets afterwards
  output="wcsim_run"
  rootname="wcsim_run"
  wcsimdir=os.environ.get('WCSIMHOME')
  particle="mu-"
  energy=1500
  nevts=1000
  split=1000
  
  try:
    myopts, args = getopt.getopt(sys.argv[1:],"o:r:p:e:n:s:w:", ["output=","rootname=","particle=","num=","energy=","split=","WCSim="])  
    print( myopts )
    print( args )
  except getopt.GetoptError, e:
    print (str(e))
    print("Usage: %s \n \t -o (--output=) name of output .mac files, without \".mac\" \n \t -r (--rootname=) name of root files to be generated, without \".root\" \n \t -p (--particle=) particle to simulate (e-, mu+ etc.) \n \t -e (--energy=) particle energy \n \t -n (--num=) total number of events to simulate \n \t -s (--split=) number of events to generate with a single mac file before splitting it \n \t -w (--WCSim=) directory containing WCSim" % sys.argv[0])
    sys.exit(2)
  for o, a in myopts:
    print( "o = " + o + "    a = " + a ) 
    if o in ("-o","--output"):
      output=a
      print( output )
    if o in ("-r","--rootname"):
      rootname=a
    if o in ("-p","--particle"):
      particle=a
    if o in ("-e","--energy"):
      try:
        energy = int(a)
      except ValueError:
        print("Energy must be an int")
        return False
    if o in ("-n","--num"):
      try:
        nevts = int(a)
      except ValueError:
        print("Number of events must be an int")
        return False
    if o in ("-s","--split"):
      try:
        split = int(a)
      except ValueError:
        print("Number of events at which to split files must be an it")
        return False
    if o in ("-w","--WCSim="):
      wcsimdir = a

  assert(split > 0)
  assert(nevts > 0)
  nFiles = int(ceil( nevts / split ))

  for file in range (0, nFiles):
    name = rootname
    name += ( "_" + str(file) )
    seed1 = file * 10
    seed2 = file * 100
    saveName = output
    saveName += ( "_" + str(file) )
    
    script = MakeScript( particle, energy, name, seed1, seed2, split )
    SaveFile( script, saveName + ".mac" )
      
    shellScript = MakeShellScript( wcsimdir, saveName + ".mac" )
    SaveFile( shellScript, saveName + ".sh" )
  
  return True


def MakeScript(particle, energy, name, seed1, seed2, num):
  myString = """
  # Sample setup macro with no visualization
  
  /run/verbose 0
  /tracking/verbose 0
  /hits/verbose 0
  
  ## select the geometry
  /WCSim/WCgeom GiantPhotonTest
  
  # uncomment to use Leigh's pmt simulation
  #/WCSim/PMTSim chips
  
  #Added for the PMT QE option 08/17/10 (XQ)
  # 1. Stacking only mean when the photon is generated
  # the QE is applied to reduce the total number of photons
  # 2. Stacking and sensitivity detector
  # In the stacking part, the maximum QE is applied to reduce 
  # the total number of photons
  # On the detector side, the rest of QE are applied according to QE/QE_max
  # distribution. This option is in particular important for the WLS
  # 3. The last option means all the QE are applied at the detector
  # Good for the low energy running. 
  /WCSim/PMTQEMethod     Stacking_Only 
  #/WCSim/PMTQEMethod     Stacking_And_SensitiveDetector
  #/WCSim/PMTQEMethod     SensitiveDetector_Only
  
  ## Whether to use the PMT collection efficiency (as a function of
  # the angle at which the incoming optical photon strikes the PMT)
  # Options: "on" or "off"
  /WCSim/PMTCollEff on
  
  /WCSim/Construct
  
  # command to choose save or not save the pi0 info 07/03/10 (XQ)
  #/WCSim/SavePi0 true
  
  ## select the input nuance-formatted vector file
  ## you can of course use your own
  #/mygen/vecfile ../genieFluxes/nue.vec
  
  
  # Or you can use the G4 General Particle Source
  # (choose 'laser' to use the particle source instead of gun)
  /mygen/generator laser
  
  #/gps/particle e-
  /gps/particle %s
  
  # can use some vtx position distribution
  # or just a point source
  /gps/pos/type Point
  /gps/pos/centre 0 0 0 cm 
  
  /gps/energy %s MeV
  
  # use direction distribution or fixed values
  #/gps/ang/type iso
  /gps/direction 0 0 1 
  
  /gps/time 0
  
  
  ## change the name of the output root file, default = wcsim.root
  #/WCSimIO/RootFile ~/some/absolute/path/file.root
  /WCSimIO/RootFile %s.root
  
  ## Whether to save an ntuple with all the optical photon tracks, default = false
  /WCSimIO/SavePhotonNtuple true
  
  ## Set the name of the photon ntuple root file, default = localfile_photons.root
  /WCSimIO/PhotonNtuple %s_photons.root
  
  # seed the random generator with two integers
  /random/setSeeds %s %s
  
  /run/beamOn %s
  #exit
  """
  return myString % ( particle, energy, name, name, seed1, seed2, num )

def MakeShellScript( wcsimdir, macName ):
  script = """#!/bin/bash
  cd %s
  . setupWCSim.sh
  WCSim %s
  echo "Done"
  """
  return ( script % (wcsimdir, macName) )


def SaveFile( script, saveName ):
  try:
    text_file = open(saveName, "w")
    text_file.write(script)
  except IOError:
    print("Error, could not open " + saveName + " for writing")
    return
  else:
    print( "Successfully wrote to " + saveName )
    text_file.close()



if __name__ == "__main__":
  main()
  
