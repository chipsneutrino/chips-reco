import sys, os, argparse
from os import listdir
import cmd

#input = "sk1pe/"
input = "pmtSim/"
#input = "tot/"

for num in range(0,10):
    cmd = ("hadd " + input + "temp" + str(num) + ".root " + input + "digitizerLikelihood_" + str(num) + "*")
    print cmd
    os.system(cmd)

cmd = ("hadd " + input + "preNorm.root " + input + "temp*")
print cmd
os.system(cmd)

#cmd = ("rm " + input + "temp* " + input + "digitizerLikelihood_*")
#print cmd
#os.system(cmd)

  


