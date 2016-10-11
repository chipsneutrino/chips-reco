=== Scattering Table Documentation ===

The scattering table aims to parametrise the amount of scattered light as a 
function of position. It is hence the only part of the reconstruction software
that is dependent on the geometry.

The simulation is used to produce a large sample of isotropic 3MeV electrons.
This is done such as to approximate a Cherenkov light point source. In fact, 
two of these samples are produced with one difference:
- One sample is simulated in the normal way.
- The second sample, refered to as "direct" contains only direct light. This
  means that scattering must be disabled.
- The python script makeScatteringJobs produces the needed files to run
  these jobs on a batch farm (will need personalising to each system).

At this point you will have a series of normal and direct files. The next step
is to produce some simplified trees from these containing the information
required by the scattering table.
- The macro runMakeNewScatteringFiles.C compiles and runs its namesake macro
  makeNewScatteringFiles.C. It takes two arguments:

root -l runMakeNewScatteringFiles.C'("/some/input/files*.root","outfile.root")'

- Wildcards accepted (and encouraged!) for the list of input files. 
- The macro needs to be run separately on the normal and direct input files.

At this point, you now have two sets of scattering trees. The next step is to 
turn these trees into a scattering table. The table itself is actually three
tables (one for the barrel, top and bottom regions) with six dimensions. It 
is clear from the dimensionality that the number of bins rapidly blows up
if you want to have n bins per dimension.

The scattering table code runs best compiled using g++. 

g++ `root-config --cflags --glibs` buildScatteringTable.C -o buildScatteringTable

It then runs taking three arugments:

./buildScatteringTable /normal/trees/*.root /direct/trees/*.root scatteringTable.root

The scatteringTable.root is then ready for use in the reconstruction software. 
Depending on statistics there might be some poorly populated bins that can 
fluctuate to having negative values. The reconstruction just assigns these to 
have some default value. It isn't really a problem, but there is an executable 
that you can use to check the table over.

g++ `root-config --cflags --glibs` checkScatteringTable.C -o checkScatteringTable

./checkScatteringTable scatteringTable.root



