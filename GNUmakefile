## Makefile ##

#
# notes:
#  WCSIM_INCDIR  points to WCSim .hh files
#  WCSIM_LIBDIR  points to WCSim .o files
#  Have added -lSpectrum to LIBS so this differs from the version WCSimAnalysis ships with -AJP 21/03/13
#

CXX						= g++
CXXDEPEND     = -MM
CXXFLAGS      = -g -O3 -Wall -fPIC
LD            = g++
LDFLAGS       = -g -O3

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLDFLAGS  := $(shell root-config --ldflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTLIBS		 += -lSpectrum

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

INCDIR = ./include
SRCDIR = ./src
TMPDIR = ./tmp
LIBDIR = ./lib
BINDIR = ./bin

INCLUDES = -I$(INCDIR)

WCSIM_INCDIR = ${WCSIMHOME}/include
WCSIM_LIBDIR = ${HOME}/geant4/tmp/Linux-g++/WCSim

WCSIM_INCLUDES = -I$(WCSIM_INCDIR)
WCSIM_LDFLAGS += -L$(WCSIM_LIBDIR)
WCSIM_LIBS += -lWCSim

.PHONY: 
	all

all: rootcint shared libWCSimAnalysis.a

ROOTSO := $(LIBDIR)/libWCSimAnalysis.so

ROOTDICT := $(SRCDIR)/WCSimAnalysisRootDict.cc

# This is the list of all the ROOT-based classes we need to worry about.
# Assumes that each class has src/*.cc, include/*.hh and tmp/*.o files.
ROOTCLASS := WCSimNvidiaMath WCSimTransmissionFunctionLookup WCSimPiZeroSeed WCSimPiZeroSeedGenerator WCSimPiZeroSeeder WCSimPiZeroHoughSeeder WCSimPiZeroSingleElectronSeeder WCSimPiZeroElectronAdjuster WCSimPiZeroFitter WCSimDetectorParameters WCSimTimePredictor WCSimTimeLikelihood2 WCSimLikelihoodTrackFactory WCSimLikelihoodTrackBase WCSimLikelihoodPhotonTrack WCSimFitterTrackParMap WCSimIntegralLookupMaker3D WCSimIntegralLookup3D WCSimIntegralLookupReader WCSimIntegralLookup WCSimIntegralLookupMaker WCSimRecoEvDisplay WCSimRecoSummary WCSimEmissionProfileManager WCSimEmissionProfiles WCSimFitterConfig WCSimFitterInterface WCSimFitterParameters WCSimFitterTree WCSimFitterPlots WCSimTimeLikelihood WCSimAnalysisConfig WCSimDigitizerLikelihood WCSimTotalLikelihood WCSimLikelihoodTrack WCSimLikelihoodDigit WCSimLikelihoodDigitArray WCSimLikelihoodTuner WCSimLikelihoodFitter WCSimDisplayViewer WCSimDisplayFactory WCSimDisplay WCSimDisplayAB WCSimEveDisplay WCSimEventWriter WCSimGeometry WCSimInterface WCSimParameters WCSimRecoObjectTable WCSimRecoFactory WCSimReco WCSimRecoAB WCSimRecoDigit WCSimRecoCluster WCSimRecoClusterDigit WCSimRecoRing WCSimRecoVertex WCSimRecoEvent WCSimTrueEvent WCSimTrueTrack WCSimHoughTransform WCSimHoughTransformArray WCSimDataCleaner WCSimVertexFinder WCSimVertexGeometry WCSimVertexViewer WCSimRingFinder WCSimRingViewer WCSimNtupleFactory WCSimNtuple WCSimRecoNtuple WCSimVertexNtuple WCSimVertexSeedNtuple WCSimNtupleWriter WCSimMsg WCSimChargePredictor WCSimRecoSeed WCSimRecoSlicer 

# Create the ROOTINC list from the class list, remembering to also add the LinkDef file
ROOTINC = $(ROOTCLASS:%=$(INCDIR)/%.hh)
ROOTINC += $(INCDIR)/WCSimAnalysisRootLinkDef.hh
# For dictionary generation, want a version without the include directory.
ROOTINCNODIR = $(ROOTCLASS:%=%.hh)
ROOTINCNODIR += WCSimAnalysisRootLinkDef.hh
# Now for the ROOTSRC list
ROOTSRC = $(ROOTCLASS:%=$(SRCDIR)/%.cc)
# Finally, the ROOTOBJ list, remembering to add on the RootDict file
ROOTOBJS = $(ROOTCLASS:%=$(TMPDIR)/%.o)
ROOTOBJS += $(TMPDIR)/WCSimAnalysisRootDict.o

ROOTEXTOBJS := $(WCSIM_LIBDIR)/WCSimRootEvent.o $(WCSIM_LIBDIR)/WCSimRootGeom.o $(WCSIM_LIBDIR)/WCSimRootDict.o

$(TMPDIR)/%.o : $(SRCDIR)/%.cc
	@echo "<**Compiling $@**>"
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) -c $< -o $@ 

$(TMPDIR)/%.d: $(SRCDIR)/%.cc
	@echo "<**Depend $@**>"
	@mkdir -p $(@D)
	set -e; $(CXX) $(CXXDEPEND) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) $< \
	| sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(ROOTDICT) : $(ROOTSRC) $(ROOTINC)
	rootcint -f $(ROOTDICT) -c -I$(INCDIR) -I$(WCSIM_INCDIR) -I$(shell root-config --incdir) $(ROOTINCNODIR)

rootcint: $(ROOTDICT)

shared: $(ROOTDICT) $(ROOTSRC) $(ROOTINC) $(ROOTOBJS)
	@echo "<**Shared**>"
	@mkdir -p lib
	$(CXX) -shared $(ROOTLIBS) $(ROOTGLIBS) -O $(ROOTOBJS) -o $(ROOTSO)

libWCSimAnalysis.a : $(ROOTOBJS)
	$(RM) $@
	ar clq $@ $(ROOTOBJS)
	
evDisp:	
	$(CXX) `root-config --cflags --glibs --libs --evelibs` -I./include ${WCSIM_INCLUDES} -L${WCSIMANAHOME} -L${WCSIMHOME} -o evDisplay evDisplay.cc ${WCSIMHOME}/src/WCSimRootDict.cc ${WCSIMANAHOME}/src/WCSimAnalysisRootDict.cc -lWCSim -lWCSimAnalysis -lEG -lSpectrum -lMinuit

fitterProfile:
	$(CXX) $(CXXFLAGS) `root-config --cflags --glibs --libs --evelibs` -I./include ${WCSIM_INCLUDES} -L${WCSIMANAHOME} -L${WCSIMHOME} -L/home/ajperch/software/gperftools/lib/ -o fitterProfile runFitterProfile.cc ${WCSIMHOME}/src/WCSimRootDict.cc ${WCSIMANAHOME}/src/WCSimAnalysisRootDict.cc -lWCSim -lWCSimAnalysis -lEG -lSpectrum -lMinuit -lprofiler

clean :
	@echo "<**Clean**>"
	rm -f $(SRCDIR)/*~ $(INCDIR)/*~ $(TMPDIR)/*.o $(TMPDIR)/*.d $(TMPDIR)/*.a $(LIBDIR)/*.so $(SRCDIR)/WCSimAnalysisRootDict.*

DEPS = $(ROOTOBJS:$(TMPDIR)/%.o=$(TMPDIR)/%.d)

ifeq ($(MAKECMDGOALS),all)
 include $(DEPS)
endif

ifeq ($(MAKECMDGOALS),shared)
 include $(DEPS)
endif
