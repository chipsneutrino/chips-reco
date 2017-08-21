#!/usr/bin/python
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#  J. Tingey - UCL (j.tingey.16@ucl.ac.uk)                                           #
#      modified from original by S. Germani - UCL  (s.germani@ucl.ac.uk)             #
#                                                                                    #
#  Script to run PID TMVA Neural Network for PID using WCSimAnalysis Rconstruction   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

import sys, os, argparse
from os import listdir
import cmd

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  parse_args     :::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def parse_args():
    parser = argparse.ArgumentParser(description='Python script to handle the running of the various stages of the PID') 
    parser.add_argument('inputDir', help = 'Path to the input directory, containing e.g. nue_cc_flux/NuMI/ElectronLike etc...')
    parser.add_argument('-o', '--outputDir', help='Path to make PID directory within, if not specified will be same as inputDir')
    parser.add_argument('-g', '--geoName', help = 'Name of the CHIPS geometry. Used for naming files.')
    parser.add_argument('-n', '--maxNum', help = 'Maximum number of files that will be used from each input dir e.g. nue_cc_flux/NuMI/ElectronLike', default = 1000)
    
    parser.add_argument('--lists', action = 'store_true') # Make the file lists that are passed to the combineTrees macro.
    parser.add_argument('--combElMu', action = 'store_true') # Combine The events for training ANNelMu
    parser.add_argument('--trainElMu', action = 'store_true') # Train the ANNelMu
    parser.add_argument('--combElNc', action = 'store_true') # Combine the events and read the ANNelMu value for the events to train the ANNelNC
    parser.add_argument('--trainElNc', action = 'store_true') # Train the ANNelNC
    parser.add_argument('--read', action = 'store_true') # Read the nue_all and numu_all events using the two trained networks
    parser.add_argument('--all', action = 'store_true') # Do all of the above...
    
    return parser.parse_args()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  checkDirsTrain     :::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def checkDirsTrain(inputDir, outputDir):
    # Check that the input directory contains all that is needed.
    print "Checking Train Directories..."      
    if not os.path.isdir(inputDir):
        print "The given input directory does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "nue_cc_flux/NuMI/ElectronLike")):
        print "The input directory 'nue_cc_flux/ElectronLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "nue_cc_flux/NuMI/MuonLike")):
        print "The input directory 'nue_cc_flux/MuonLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "numu_cc_flux/NuMI/ElectronLike")):
        print "The input directory 'numu_cc_flux/ElectronLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "numu_cc_flux/NuMI/MuonLike")):
        print "The input directory 'numu_cc_flux/MuonLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "numu_nc_flux/NuMI/ElectronLike")):
        print "The input directory 'numu_nc_flux/ElectronLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "numu_nc_flux/NuMI/MuonLike")):
        print "The input directory 'numu_nc_flux/MuonLike' does not exist: Exit!"
        sys.exit(1)
              
    # Now create the directory to store the PID files...
    if not os.path.isdir(os.path.join(outputDir, "PID/")):
        os.mkdir(os.path.join(outputDir, "PID/"))            
        
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  checkDirsRead     ::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def checkDirsRead(inputDir, outputDir):
    # Check that the input directory contains all that is needed.
    print "Checking Read Directories..."      
    if not os.path.isdir(os.path.join(inputDir, "nue_all_flux/NuMI/ElectronLike")):
        print "The input directory 'nue_all_flux/ElectronLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "nue_all_flux/NuMI/MuonLike")):
        print "The input directory 'nue_all_flux/MuonLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "numu_all_flux/NuMI/ElectronLike")):
        print "The input directory 'numu_all_flux/ElectronLike' does not exist: Exit!"
        sys.exit(1)
    if not os.path.isdir(os.path.join(inputDir, "numu_all_flux/NuMI/MuonLike")):
        print "The input directory 'numu_all_flux/MuonLike' does not exist: Exit!"
        sys.exit(1)
        
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  makeLists   ::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
def makeLists(inputDir, outputDir, geoName, type, maxNum):
    print "Making lists for type -> " + type + "..."
    if not os.path.isdir(os.path.join(outputDir, "PID/fileLists")):
        os.mkdir(os.path.join(outputDir, "PID/fileLists"))
    # check that the directories have approx same number of files and decide how to split for training/sample 
    elDir = os.path.join(inputDir, type + "/NuMI/ElectronLike")
    muDir = os.path.join(inputDir, type + "/NuMI/MuonLike")
    num_el_files = len([f for f in os.listdir(elDir)
                        if f.endswith("tree.root")])
    num_mu_files = len([f for f in os.listdir(muDir)
                        if f.endswith("tree.root")])
    print "Num files ElectronLike -> " + str(num_el_files) + ", Num files MuonLike -> " + str(num_mu_files)
    if (num_el_files / num_mu_files) > 0.95 or (num_el_files / num_mu_files) < 1.05:
        print 'Accepted. Number of files are similar...'
    else:
        print 'Rejected. Number of files differ to much...'
        exit(1)
        
    if type == "nue_cc_flux":
        actualType = "nuel_cc"
    elif type == "numu_cc_flux":
        actualType = "numu_cc"
    elif type == "numu_nc_flux":
        actualType = "numu_nc"
    elif type == "nue_all_flux":
        actualType = "nuelAll"
    elif type == "numu_all_flux":
        actualType = "numuAll"
        
    # Open Files for both training and sample file lists...
    el_File = open(os.path.join(outputDir, "PID/fileLists/" + geoName + "_" + actualType + "_elFit_List.txt"), "w")
    mu_File = open(os.path.join(outputDir, "PID/fileLists/" + geoName + "_" + actualType + "_muFit_List.txt"), "w")
    sampleCount = 0
    #Loop through dirs and match files together etc...     
    for elF in listdir(os.path.join(inputDir, type + "/NuMI/ElectronLike")):
        for muF in listdir(os.path.join(inputDir, type + "/NuMI/MuonLike")):
            if elF.endswith("tree.root") and muF.endswith("tree.root"):
                if elF[-33:-30] == muF[-33:-30] and elF[-29:-25] == muF[-29:-25]: #The new one
                    el_File.write(os.path.join(os.path.join(inputDir, type + "/NuMI/ElectronLike"), elF) + "\n")
                    mu_File.write(os.path.join(os.path.join(inputDir, type + "/NuMI/MuonLike"), muF) + "\n") 
                    sampleCount += 1
                        
            if int(sampleCount) >= int(maxNum):
                break      
        if int(sampleCount) >= int(maxNum):
            break                             
                    
    el_File.close()   
    mu_File.close() 
    print "Files Used -> " + str(sampleCount)
    
    
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  wc_CombineTrees     ::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
def run_wc_CombineTrees(outputDir, type, geoName):
    fileListDir = os.path.join(outputDir, "PID/fileLists")
    if not os.path.isdir(os.path.join(outputDir, "PID/combinedFiles")):
        os.mkdir(os.path.join(outputDir, "PID/combinedFiles"))
    #Find the appropriate files
    for file in listdir(fileListDir):
        if file.endswith(".txt"):
            if file[-22:-15] == type and file[-14:-9] == "muFit":
                muFile = os.path.join(fileListDir, file);
            if file[-22:-15] == type and file[-14:-9] == "elFit":
                elFile = os.path.join(fileListDir, file);

    tempName = os.path.join(outputDir, "PID/combinedFiles/") + "temp_PIDTree_"   
    #Set root command and run it
    cmd = ("root -l -q \'$WCSIMANAHOME/pid/macros/wc_combineTrees.C(\"%s\",\"%s\", \"%s\")\'" % (muFile, elFile, tempName))
    os.system(cmd)
        
    counter = 0   
    for hadFile in listdir(os.path.join(outputDir, "PID/combinedFiles/")):
        if hadFile[-8:] == "ann.root":
            counter += 1
            
    print "Num Files to be Hadded-> " + str(counter) 
    
    if counter > 1:
        cmd = ("hadd -f " + os.path.join(outputDir, "PID/combinedFiles/") + geoName + "_" + type + "_combined.root" + " " + os.path.join(outputDir, "PID/combinedFiles/") + "*ann.root")
        os.system(cmd)      
        cmd = ("rm " + os.path.join(outputDir, "PID/combinedFiles/") + "*ann.root")
        os.system(cmd)     
    else:
        cmd = ("mv " + os.path.join(outputDir, "PID/combinedFiles/") + "*ann.root" + " " + os.path.join(outputDir, "PID/combinedFiles/") + geoName + "_" + type + "_combined.root")
        os.system(cmd)
           
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  run_wctrainTMVA_electronMuonCCQE  ::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
def run_wc_trainTMVA_electronMuonCCQE(outputDir):
    if not os.path.isdir(os.path.join(outputDir, "PID/weights")):
        os.mkdir(os.path.join(outputDir, "PID/weights"))
    combFilesDir = os.path.join(outputDir, "PID/combinedFiles/")
    weightsDir = os.path.join(outputDir, "PID/")
    for file in listdir(combFilesDir):
        if file.endswith("combined.root"):
            if file[-21:-14] == "nuel_cc":
                el_cc_file = file
            if file[-21:-14] == "numu_cc":
                mu_cc_file = file
                
    cmd = ("root -l -q \'$WCSIMANAHOME/pid/macros/wc_trainTMVA_electronMuonCCQE.C(\"%s\",\"%s\",\"%s\")\'" % (os.path.join(combFilesDir, el_cc_file), os.path.join(combFilesDir, mu_cc_file), weightsDir)) 
    os.system(cmd) 
    
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  run_wc_readTMVA_NCSelect     :::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
def run_wc_readTMVA_NCSelect(outputDir, type, geoName):
    fileListDir = os.path.join(outputDir, "PID/fileLists")
    pidDir = os.path.join(outputDir, "PID")
    #Find the appropriate files
    for file in listdir(fileListDir):
        if file.endswith(".txt"):
            if file[-22:-15] == type and file[-14:-9] == "muFit":
                muFile = os.path.join(fileListDir, file);
            if file[-22:-15] == type and file[-14:-9] == "elFit":
                elFile = os.path.join(fileListDir, file);

    tempName = os.path.join(outputDir, "PID/combinedFiles/") + "temp_PIDTree_"   
    #Set root command and run it
    cmd = ("root -l -q \'$WCSIMANAHOME/pid/macros/wc_readTMVA_NCSelect.C(\"%s\",\"%s\",\"%s\",\"%s\")\'" % (muFile, elFile, tempName, pidDir))
    os.system(cmd)
        
    counter = 0   
    for hadFile in listdir(os.path.join(outputDir, "PID/combinedFiles/")):
        if hadFile[-8:] == "ann.root":
            counter += 1
            
    print "Num Files to be Hadded-> " + str(counter) 
    
    if counter > 1:
        cmd = ("hadd -f " + os.path.join(outputDir, "PID/combinedFiles/") + geoName + "_" + type + "_combined_NCSelect.root" + " " + os.path.join(outputDir, "PID/combinedFiles/") + "*ann.root")
        os.system(cmd)      
        cmd = ("rm " + os.path.join(outputDir, "PID/combinedFiles/") + "*ann.root")
        os.system(cmd)     
    else:
        cmd = ("mv " + os.path.join(outputDir, "PID/combinedFiles/") + "*ann.root" + " " + os.path.join(outputDir, "PID/combinedFiles/") + geoName + "_" + type + "_combined_NCSelect.root")
        os.system(cmd)           

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  run_wcTrainTMVA  :::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
def run_wc_trainTMVA_electronCCQEvsNC(outputDir):
    if not os.path.isdir(os.path.join(outputDir, "PID/weights")):
        os.mkdir(os.path.join(outputDir, "PID/weights"))
    combFilesDir = os.path.join(outputDir, "PID/combinedFiles/")
    weightsDir = os.path.join(outputDir, "PID/")
    for file in listdir(combFilesDir):
        if file.endswith("combined_NCSelect.root"):
            if file[-30:-23] == "nuel_cc":
                el_cc_file = file
            if file[-30:-23] == "numu_nc":    
                mu_nc_file = file
                         
    cmd = ("root -l -q \'$WCSIMANAHOME/pid/macros/wc_trainTMVA_electronCCQEvsNC.C(\"%s\",\"%s\",\"%s\")\'" % (os.path.join(combFilesDir, el_cc_file), os.path.join(combFilesDir, mu_nc_file), weightsDir))
    os.system(cmd)   
    
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  wc_readTMVA       ::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
def run_wcReadTMVA(outputDir, type, geoName):
    if not os.path.isdir(os.path.join(outputDir, "PID/output")):
        os.mkdir(os.path.join(outputDir, "PID/output"))
    fileListDir = os.path.join(outputDir, "PID/fileLists")
    pidDir = os.path.join(outputDir, "PID")
    #Find the appropriate files
    for file in listdir(fileListDir):
        if file.endswith(".txt"):
            if file[-22:-15] == type and file[-14:-9] == "muFit":
                muFile = os.path.join(fileListDir, file);
            if file[-22:-15] == type and file[-14:-9] == "elFit":
                elFile = os.path.join(fileListDir, file);

    tempName = os.path.join(outputDir, "PID/output/") + "temp_PIDTree_"   
    #Set root command and run it
    cmd = ("root -l -q \'$WCSIMANAHOME/pid/macros/wc_readTMVA.C(\"%s\",\"%s\",\"%s\",\"%s\")\'" % (muFile, elFile, tempName, pidDir))
    os.system(cmd)
        
    counter = 0   
    for hadFile in listdir(os.path.join(outputDir, "PID/output/")):
        if hadFile[-8:] == "ann.root":
            counter += 1
            
    print "Num Files to be Hadded-> " + str(counter) 
    
    if counter > 1:
        cmd = ("hadd -f " + os.path.join(outputDir, "PID/output/") + geoName + "_" + type + "_readOutput.root" + " " + os.path.join(outputDir, "PID/output/") + "*ann.root")
        os.system(cmd)      
        cmd = ("rm " + os.path.join(outputDir, "PID/output/") + "*ann.root")
        os.system(cmd)     
    else:
        cmd = ("mv " + os.path.join(outputDir, "PID/output/") + "*ann.root" + " " + os.path.join(outputDir, "PID/output/") + geoName + "_" + type + "_readOutput.root")
        os.system(cmd)    
    
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  run_wcPID      :::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def run_wcPID():
    
    print "#### RUN THE PID ####"
    
    # Check Arguments
    args = parse_args()
    if not args:
        print 'Invalid Arguments'
        exit(1)
        
    if args.outputDir:
        outputDir = args.outputDir
    else:
        outputDir = args.inputDir
        
    print "Input Path -> " + args.inputDir
    print "Output path-> " + outputDir
    
    # Check the directories needed to train the PID
    checkDirsTrain(args.inputDir, outputDir)
    
    if args.lists or args.all:
        print "#### MAKING THE TRAINING FILE LISTS ####"
        makeLists(args.inputDir, outputDir, args.geoName, "nue_cc_flux", args.maxNum)
        makeLists(args.inputDir, outputDir, args.geoName, "numu_cc_flux", args.maxNum)
        makeLists(args.inputDir, outputDir, args.geoName, "numu_nc_flux", args.maxNum)
    
    if args.combElMu or args.all:
        print "#### COMBINING ELECTRON/MUON FILES ####"
        run_wc_CombineTrees(outputDir, "nuel_cc", args.geoName)
        run_wc_CombineTrees(outputDir, "numu_cc", args.geoName)
        
    if args.trainElMu or args.all:
        print "#### TRAINING THE ELECTRON/MUON NEURAL NETWORK ####"
        run_wc_trainTMVA_electronMuonCCQE(outputDir)
        
    if args.combElNc or args.all:
        print "#### COMBINING ELECTRON/NC FILES + READING ANNElMu VALUE ####" 
        run_wc_readTMVA_NCSelect(outputDir, "nuel_cc", args.geoName)
        run_wc_readTMVA_NCSelect(outputDir, "numu_nc", args.geoName)
    
    if args.trainElNc or args.all:
        print "#### TRAINING THE ELECTRON/NC NEURAL NETWORK ####"
        run_wc_trainTMVA_electronCCQEvsNC(outputDir)
        
    if args.read or args.all:
        checkDirsRead(args.inputDir, outputDir)
        print "#### MAKING THE SAMPLE FILE LISTS ####"
        makeLists(args.inputDir, outputDir, args.geoName, "nue_all_flux", args.maxNum)
        makeLists(args.inputDir, outputDir, args.geoName, "numu_all_flux", args.maxNum)
        print "#### READING THE SAMPLE EVENTS ####"
        run_wcReadTMVA(outputDir, "nuelAll", args.geoName)
        run_wcReadTMVA(outputDir, "numuAll", args.geoName)
        
    print "#### FINISHED RUNNING THE PID ####" 
        
if __name__ == "__main__":
   run_wcPID()
