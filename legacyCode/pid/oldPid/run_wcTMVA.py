#!/usr/bin/python
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#  S. Germani - UCL  (s.germani@ucl.ac.uk)
#
#  Script to run PID TMVA Neural Network for  PID using  WCSimAnalysis Rconstruction 
#  Arguments are:
#       1: Muon     Hypothesis Files List
#       2: Electron Hypothesis Files List
#       3: Output Directory (Optional - Default is ./)
#       4: Common String to add to the Output Files (Optional)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

import sys,os

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  run_wcTMVA      ::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def run_wcTMVA(argv):
    """Run TMVA Neural Network for a list of Muon and Electron hypotheses files from  WCSimAnalysis Rconstruction
       arguments are:
       1: Muon     Hypothesis Files List
       2: Electron Hypothesis Files List
       3: Output Directory (Optional - Default is ./)
       4: Common String to add to the Output Files (Optional)"""


    #:::::::::::::::  Check Arguments   :::::::::::::::::::::::::::::::::::::::::::
    if len(argv) < 2:
        print "usage: run_wcTMVA.py   mu_list  el_list   [output_dir] [common_name]"
        exit(1)

    for arg in argv:
        print arg

    #  Get Inout Files Lists and Output dir   :
    mu_list_file = argv[0]
    el_list_file = argv[1]
    out_dir     = "./"
    common_name = ""

    if len(argv)>=3 :
        out_dir = argv[2]
        
    if len(argv)>=4 :
        common_name = argv[3]

    #:::::::::::::::  Open Files Lists   :::::::::::::::::::::::::::::::::::::::::::
    mu_list = open(mu_list_file, 'r')
    el_list = open(el_list_file, 'r')
    
    mu_lines = mu_list.readlines()
    el_lines = el_list.readlines()
    
    
    #:::::::::::::::  Loop over Files Lists   ::::::::::::::::::::::::::::::::::::::
    for mu_line, el_line  in zip(mu_lines, el_lines):
        muname = os.path.basename(mu_line.split("\n")[0])
        elname = os.path.basename(el_line.split("\n")[0])
        prefix = os.path.commonprefix([muname, elname])
        
        #Set output file name
        outname = out_dir+"/"+common_name+prefix+"PIDTree_ann.root"
        
        #Set root command and run it
        cmd = ("root -l -q \'$WCSIMANAHOME/pid/run_wcTMVA.C(\"%s\",\"%s\", \"%s\")\'" % (mu_line.split("\n")[0], el_line.split("\n")[0], outname) )
        print cmd
        os.system(cmd)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::  Main            ::::::::::::::::::::::::::::::::::::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if __name__ == "__main__":
   run_wcTMVA(sys.argv[1:])
