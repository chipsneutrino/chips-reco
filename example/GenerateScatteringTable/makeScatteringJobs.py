#! /usr/bin/env python

# Python script to generate isotropic 3MeV electrons
import math
import random

# ============================================================================
# Generate random position within a circle of radius r
def GeneratePositionCircle(r):

  position = [0,0]

  u = random.random()*2*math.pi
  v = math.sqrt(random.random())

  position[0] = r * v * math.cos(u)
  position[1] = r * v * math.sin(u)

  return position

# ============================================================================
# Generate the random isotropic direction
def GenerateIsotropicDirection():
  direction = [0,0,0]

  # Get two random numbers between 0 and 1
  u = random.random()
  v = random.random()

  # Use these to get theta and phi
  theta = 2*math.pi*u
  phi = math.acos((2*v)-1)
  
  # Now convert these into x, y and z
  direction[0] = math.sin(theta)*math.cos(phi)
  direction[1] = math.sin(theta)*math.sin(phi)
  direction[2] = math.cos(theta)
  
  # We have our direction, return it!
  return direction

# ============================================================================

# Main script

# Electron energy in MeV
energy = 3
# Number of jobs
nvec = 800
# Number of events per vector file
nevents = 10000
# Number of particles per bomb
# nelec = 1
# Particle pdg code
pdgCode = 11
# Detector radius in cm
dRadius = 600 # Half of 1200
# Detector half height in cm
dHeight = 485 # Half of 970

# For each vec file that we produce, we also need the following:
# - A shell script to launch the job with scattering.
# - A shell script to launch the job with direct light only.
# - A mac file to run the named job in WCSim.
# Globally speaking, we then need two other sets of files:
# - One geoSetup.mac containing the geometry that we wish to use.
# - One jobOptions.mac and jobOptions2.mac
# - One version of tuning_parameters.mac with scattering on.
# - One version of tuning_parameters.mac with scattering off.

# Global stuff first
# Tuning parameters file for normal running
tunenormname = "/unix/fnu/lwhitehead/scatteringScripts/mac/tuning_parameters_normal.mac"
file = open(tunenormname,"w")
file.write("/WCSim/tuning/abwff 100000.0\n")
file.write("/WCSim/tuning/rgcff 0.24\n")
file.write("/WCSim/tuning/rayff 0.8\n")
file.write("/WCSim/tuning/bsrff 0.7\n")
file.write("/WCSim/tuning/mieff 0.0\n")
file.write("/WCSim/tuning/topveto 0\n")
file.write("/WCSim/tuning/tvspacing 100\n")
file.close()

# Tuning parametesr for direct light only
tunedirectname = "/unix/fnu/lwhitehead/scatteringScripts/mac/tuning_parameters_direct.mac"
file = open(tunedirectname,"w")
file.write("/WCSim/tuning/abwff 100000.0\n")
file.write("/WCSim/tuning/rgcff 0.24\n")
file.write("/WCSim/tuning/rayff 100000.0\n")
file.write("/WCSim/tuning/bsrff 0.0\n")
file.write("/WCSim/tuning/mieff 100000.0\n")
file.write("/WCSim/tuning/topveto 0\n")
file.write("/WCSim/tuning/tvspacing 100\n")
file.close()

# Geosetup mac file
geoname = "/unix/fnu/lwhitehead/scatteringScripts/mac/geoSetup.mac"
file = open(geoname,"w")
file.write("/WCSim/WCgeom CHIPS_10kton_10inch_HQE_10perCent\n")
file.write("/WCSim/PMTQEMethod Stacking_And_SensitiveDetector\n")
file.close()

jobname = "/unix/fnu/lwhitehead/scatteringScripts/mac/jobOptions.mac"
file = open(jobname,"w")
file.write("/WCSim/physics/list WCSim\n")
file.close()

job2name = "/unix/fnu/lwhitehead/scatteringScripts/mac/jobOptions2.mac"
file = open(job2name,"w")
file.write("/WCSim/physics/secondaries/model BINARY\n")
file.close()

submitname = "/home/lwhitehead/work/CHIPS/scattering/submitJobs.sh"
jobfile = open(submitname,"w")
jobfile.write("#! /bin/bash/\n")

for vecfiles in range (0,nvec):

  # Set the random seed
  random.seed(vecfiles)
  
  # Create the vector file
  particleName = ""
  mass=0.
  if(pdgCode == 11):
    particleName = "Electron"
    mass = 0.511
  else:
    particleName = "Muon"
    mass = 105.658

#  vecfilename = "/unix/fnu/lwhitehead/scatteringScripts/electronBomb_" + str(vecfiles) + ".vec"
#  file = open(vecfilename,"w")
#  
#  for events in range(0,nevents):
#    file.write("$ begin\n")
#    file.write("$ nuance "+str(events)+"\n")
#    
#    # Generate a random position inside a cylinder (x,y)
#    tempPos = GeneratePositionCircle(dRadius)
#    # Now generate the random z position in the range -dHeight to dHeight
#    tempZ = random.random() * (2*dHeight) - dHeight;
#    file.write("$ vertex " + str(tempPos[0]) + " " + str(tempPos[1]) + " " + str(tempZ) + " 0\n")
#    # Let's make a fake neutrino and target
#    file.write("$ track 12 5.0 0.0 0.0 1.0 -1\n")
#    file.write("$ track 8016 3.0 -999 -999 -999 -1\n")
#    
#    for i in range(0,nelec):
#      tempDir = GenerateIsotropicDirection()
#      file.write("$ track " + str(pdgCode) + " " +str(energy+mass) + " " + str(tempDir[0]) + " " + str(tempDir[1]) + " " + str(tempDir[2]) + " 0\n" )
#
#		# End the event  
#    file.write("$ end\n")
#
#	# Close the current file and report success
#  file.close()
#
#  if (vecfiles % 100 == 0):
#    print "File " + vecfilename + " written."

  # Now that we have the vector file, make the files to go with it.

  # This is the normal version with all light on.
  macnormname = "/unix/fnu/lwhitehead/scatteringScripts/scatteringTable_normal_" + str(vecfiles) + ".mac"
  outnormname = "/unix/fnu/lwhitehead/scatteringFiles/scatteringTable_normal_" + str(vecfiles) + ".root"

  file = open(macnormname,"w")
  file.write("/run/verbose 0\n")
  file.write("/tracking/verbose 0\n")
  file.write("/hits/verbose 0\n")
  file.write("/mygen/generator gps\n")
  file.write("/gps/particle e-\n")
  file.write("/gps/pos/type Volume\n")
  file.write("/gps/pos/shape Cylinder\n")
  file.write("/gps/pos/radius 12.0 m\n")
  file.write("/gps/pos/halfz 9.5 m\n")
  file.write("/gps/energy 3 MeV\n")
  file.write("/gps/ang/type iso\n")
  file.write("/gps/time 0\n")
  file.write("/WCSimIO/RootFile scatteringTable_normal_" + str(vecfiles) + ".root" + "\n")
  file.write("/WCSimIO/SavePhotonNtuple false\n")
  file.write("/WCSimTrack/PercentCherenkovPhotonsToDraw 0.0\n")
  file.write("/random/setSeeds " + str(vecfiles+1) +  " " + str(2*(vecfiles+1)) + "\n")
  file.write("/run/beamOn " + str(nevents) + "\n")
  file.close()

  # This is the version with only direct light.
  macdirectname = "/unix/fnu/lwhitehead/scatteringScripts/scatteringTable_direct_" + str(vecfiles) + ".mac"
  outdirectname = "/unix/fnu/lwhitehead/scatteringFiles/scatteringTable_direct_" + str(vecfiles) + ".root"

  file = open(macdirectname,"w")
  file.write("/run/verbose 0\n")
  file.write("/tracking/verbose 0\n")
  file.write("/hits/verbose 0\n")
  file.write("/mygen/generator gps\n")
  file.write("/gps/particle e-\n")
  file.write("/gps/pos/type Volume\n")
  file.write("/gps/pos/shape Cylinder\n")
  file.write("/gps/pos/radius 12.0 m\n")
  file.write("/gps/pos/halfz 9.5 m\n")
  file.write("/gps/energy 3 MeV\n")
  file.write("/gps/ang/type iso\n")
  file.write("/gps/time 0\n")
  file.write("/WCSimIO/RootFile scatteringTable_direct_" + str(vecfiles) + ".root" + "\n")
  file.write("/WCSimIO/SavePhotonNtuple false\n")
  file.write("/WCSimTrack/PercentCherenkovPhotonsToDraw 0.0\n")
  file.write("/random/setSeeds " + str(vecfiles+1) +  " " + str(2*(vecfiles+1)) + "\n")
  file.write("/run/beamOn " + str(nevents) + "\n")
  file.close()

  # Finally, make the shell script to run the job
  jobscriptname = "/unix/fnu/lwhitehead/scatteringScripts/runScatteringJob_" + str(vecfiles) + ".sh"
  jobfile.write("qsub -q medium " + jobscriptname + "\n")
  file = open(jobscriptname,"w")
  file.write("#!/bin/bash\n")
  file.write("CURDIR=`pwd`\n")
  file.write("JOBDIR=scatterJob_" + str(vecfiles) + "\n")
  file.write("cd /tmp/ # Worker node scratch space\n")
  file.write("mkdir $JOBDIR\n")
  file.write("cd $JOBDIR\n")
  file.write(". ~/work/CHIPS/recon/test/WCSim/setupScripts/setupWCSimUCL.sh\n");
  file.write("# Copy some files across\n");
  file.write("cp " + geoname + " ./\n")
  file.write("cp " + macnormname + " normal.mac\n")
  file.write("cp " + macdirectname + " direct.mac\n")
  file.write("cp " + jobname + " ./\n")
  file.write("cp " + job2name + " ./\n")
  file.write("cp " + tunenormname + " tuning_parameters.mac\n")
  file.write("WCSim normal.mac\n")
  file.write("cp " + tunedirectname + " tuning_parameters.mac\n")
  file.write("WCSim direct.mac\n")
  file.write("mv *.root /unix/fnu/lwhitehead/scatteringFiles/\n")
  file.write("cd $CURDIR\n")
  file.write("rm -r /tmp/$JOBDIR\n")

# Now finish the shell script to launch the jobs
jobfile.close()

# The end.


