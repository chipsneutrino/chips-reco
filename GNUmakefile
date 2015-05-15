## Makefile ##

#
# notes:
#  WCSIM_INCDIR  points to WCSim .hh files
#  WCSIM_LIBDIR  points to WCSim .o files
#  Have added -lSpectrum to LIBS so this differs from the version WCSimAnalysis ships with -AJP 21/03/13
#

CXX           = g++
CXXDEPEND     = -MM
CXXFLAGS      = -g -Wall -fPIC
LD            = g++
LDFLAGS       = -g 

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

all: rootcint lib shared libWCSimAnalysis.a

ROOTSO := $(LIBDIR)/libWCSimAnalysis.so

ROOTDICT := $(SRCDIR)/WCSimAnalysisRootDict.cc

ROOTSRC := $(SRCDIR)/WCSimIntegralLookupMaker3D.cc $(INCDIR)/WCSimIntegralLookupMaker3D.hh $(SRCDIR)/WCSimIntegralLookup3D.cc $(INCDIR)/WCSimIntegralLookup3D.hh $(SRCDIR)/WCSimIntegralLookupReader.cc $(INCDIR)/WCSimIntegralLookupReader.hh $(SRCDIR)/WCSimIntegralLookup.cc $(INCDIR)/WCSimIntegralLookup.hh $(SRCDIR)/WCSimIntegralLookupMaker.cc $(INCDIR)/WCSimIntegralLookupMaker.hh $(SRCDIR)/WCSimRecoEvDisplay.cc ${INCDIR}/WCSimRecoEvDisplay.hh ${SRCDIR}/WCSimRecoSummary.cc ${INCDIR}/WCSimRecoSummary.hh ${SRCDIR}/WCSimEmissionProfiles.cc ${INCDIR}/WCSimEmissionProfiles.hh ${SRCDIR}/WCSimFitterConfig.cc ${INCDIR}/WCSimFitterConfig.hh ${SRCDIR}/WCSimFitterInterface.cc ${INCDIR}/WCSimFitterInterface.hh ${SRCDIR}/WCSimFitterParameters.cc ${INCDIR}/WCSimFitterParameters.hh ${SRCDIR}/WCSimFitterTree.cc ${INCDIR}/WCSimFitterTree.hh ${SRCDIR}/WCSimFitterPlots.cc ${INCDIR}/WCSimFitterPlots.hh ${SRCDIR}/WCSimTimeLikelihood.cc ${INCDIR}/WCSimTimeLikelihood.hh ${SRCDIR}/WCSimAnalysisConfig.cc ${INCDIR}/WCSimAnalysisConfig.hh ${SRCDIR}/WCSimDigitizerLikelihood.cc ${INCDIR}/WCSimDigitizerLikelihood.hh $(SRCDIR)/WCSimTotalLikelihood.cc $(INCDIR)/WCSimTotalLikelihood.hh $(SRCDIR)/WCSimLikelihoodTrack.cc $(INCDIR)/WCSimLikelihoodTrack.hh $(SRCDIR)/WCSimLikelihoodDigit.cc $(INCDIR)/WCSimLikelihoodDigit.hh  $(SRCDIR)/WCSimLikelihoodDigitArray.cc $(INCDIR)/WCSimLikelihoodDigitArray.hh $(SRCDIR)/WCSimLikelihoodTuner.cc $(INCDIR)/WCSimLikelihoodTuner.hh $(SRCDIR)/WCSimLikelihoodFitter.cc $(INCDIR)/WCSimLikelihoodFitter.hh $(SRCDIR)/WCSimDisplayViewer.cc $(INCDIR)/WCSimDisplayViewer.hh $(SRCDIR)/WCSimDisplayFactory.cc $(INCDIR)/WCSimDisplayFactory.hh $(SRCDIR)/WCSimDisplay.cc $(INCDIR)/WCSimDisplay.hh $(SRCDIR)/WCSimDisplayAB.cc $(INCDIR)/WCSimDisplayAB.hh $(SRCDIR)/WCSimEveDisplay.cc $(INCDIR)/WCSimEveDisplay.hh $(SRCDIR)/WCSimEventWriter.cc $(INCDIR)/WCSimEventWriter.hh $(SRCDIR)/WCSimGeometry.cc $(INCDIR)/WCSimGeometry.hh $(SRCDIR)/WCSimInterface.cc $(INCDIR)/WCSimInterface.hh $(SRCDIR)/WCSimParameters.cc $(INCDIR)/WCSimParameters.hh $(SRCDIR)/WCSimRecoObjectTable.cc $(INCDIR)/WCSimRecoObjectTable.hh $(SRCDIR)/WCSimRecoFactory.cc $(INCDIR)/WCSimRecoFactory.hh $(SRCDIR)/WCSimReco.cc $(INCDIR)/WCSimReco.hh $(SRCDIR)/WCSimRecoAB.cc $(INCDIR)/WCSimRecoAB.hh $(SRCDIR)/WCSimRecoDigit.cc $(INCDIR)/WCSimRecoDigit.hh $(SRCDIR)/WCSimRecoCluster.cc $(INCDIR)/WCSimRecoCluster.hh $(SRCDIR)/WCSimRecoClusterDigit.cc $(INCDIR)/WCSimRecoClusterDigit.hh $(SRCDIR)/WCSimRecoRing.cc $(INCDIR)/WCSimRecoRing.hh $(SRCDIR)/WCSimRecoVertex.cc $(INCDIR)/WCSimRecoVertex.hh $(SRCDIR)/WCSimRecoEvent.cc $(INCDIR)/WCSimRecoEvent.hh $(SRCDIR)/WCSimTrueEvent.cc $(INCDIR)/WCSimTrueEvent.hh $(SRCDIR)/WCSimTrueTrack.cc $(INCDIR)/WCSimTrueTrack.hh $(SRCDIR)/WCSimHoughTransform.cc $(INCDIR)/WCSimHoughTransform.hh $(SRCDIR)/WCSimHoughTransformArray.cc $(INCDIR)/WCSimHoughTransformArray.hh $(SRCDIR)/WCSimDataCleaner.cc $(INCDIR)/WCSimDataCleaner.hh $(SRCDIR)/WCSimVertexFinder.cc $(INCDIR)/WCSimVertexFinder.hh $(SRCDIR)/WCSimVertexGeometry.cc $(INCDIR)/WCSimVertexGeometry.hh $(SRCDIR)/WCSimVertexViewer.cc $(INCDIR)/WCSimVertexViewer.hh $(SRCDIR)/WCSimRingFinder.cc $(INCDIR)/WCSimRingFinder.hh $(SRCDIR)/WCSimRingViewer.cc $(INCDIR)/WCSimRingViewer.hh $(SRCDIR)/WCSimNtupleFactory.cc $(INCDIR)/WCSimNtupleFactory.hh $(SRCDIR)/WCSimNtuple.cc $(INCDIR)/WCSimNtuple.hh $(SRCDIR)/WCSimRecoNtuple.cc $(INCDIR)/WCSimRecoNtuple.hh $(SRCDIR)/WCSimVertexNtuple.cc $(INCDIR)/WCSimVertexNtuple.hh $(SRCDIR)/WCSimVertexSeedNtuple.cc $(INCDIR)/WCSimVertexSeedNtuple.hh $(SRCDIR)/WCSimNtupleWriter.cc $(INCDIR)/WCSimNtupleWriter.hh $(SRCDIR)/WCSimMsg.cc $(INCDIR)/WCSimMsg.hh $(SRCDIR)/WCSimChargePredictor.cc $(INCDIR)/WCSimChargePredictor.hh $(INCDIR)/WCSimAnalysisRootLinkDef.hh

ROOTOBJS := $(TMPDIR)/WCSimIntegralLookupMaker3D.o $(TMPDIR)/WCSimIntegralLookup3D.o $(TMPDIR)/WCSimIntegralLookupReader.o $(TMPDIR)/WCSimIntegralLookup.o $(TMPDIR)/WCSimIntegralLookupMaker.o ${TMPDIR}/WCSimRecoEvDisplay.o ${TMPDIR}/WCSimRecoSummary.o ${TMPDIR}/WCSimEmissionProfiles.o ${TMPDIR}/WCSimFitterConfig.o ${TMPDIR}/WCSimFitterInterface.o ${TMPDIR}/WCSimFitterParameters.o ${TMPDIR}/WCSimFitterTree.o ${TMPDIR}/WCSimFitterPlots.o ${TMPDIR}/WCSimTimeLikelihood.o ${TMPDIR}/WCSimAnalysisConfig.o ${TMPDIR}/WCSimDigitizerLikelihood.o $(TMPDIR)/WCSimTotalLikelihood.o $(TMPDIR)/WCSimLikelihoodTrack.o $(TMPDIR)/WCSimLikelihoodDigit.o $(TMPDIR)/WCSimLikelihoodDigitArray.o $(TMPDIR)/WCSimChargePredictor.o $(TMPDIR)/WCSimLikelihoodTuner.o $(TMPDIR)/WCSimLikelihoodFitter.o $(TMPDIR)/WCSimDisplayViewer.o $(TMPDIR)/WCSimDisplayFactory.o $(TMPDIR)/WCSimDisplay.o $(TMPDIR)/WCSimDisplayAB.o $(TMPDIR)/WCSimEveDisplay.o $(TMPDIR)/WCSimEventWriter.o $(TMPDIR)/WCSimGeometry.o $(TMPDIR)/WCSimInterface.o $(TMPDIR)/WCSimParameters.o $(TMPDIR)/WCSimRecoObjectTable.o $(TMPDIR)/WCSimRecoFactory.o $(TMPDIR)/WCSimReco.o $(TMPDIR)/WCSimRecoAB.o $(TMPDIR)/WCSimRecoDigit.o $(TMPDIR)/WCSimRecoCluster.o $(TMPDIR)/WCSimRecoClusterDigit.o $(TMPDIR)/WCSimRecoRing.o $(TMPDIR)/WCSimRecoVertex.o $(TMPDIR)/WCSimRecoEvent.o $(TMPDIR)/WCSimTrueEvent.o $(TMPDIR)/WCSimTrueTrack.o $(TMPDIR)/WCSimHoughTransform.o $(TMPDIR)/WCSimHoughTransformArray.o $(TMPDIR)/WCSimDataCleaner.o $(TMPDIR)/WCSimVertexFinder.o $(TMPDIR)/WCSimVertexGeometry.o $(TMPDIR)/WCSimVertexViewer.o $(TMPDIR)/WCSimRingFinder.o $(TMPDIR)/WCSimRingViewer.o $(TMPDIR)/WCSimNtupleFactory.o $(TMPDIR)/WCSimNtuple.o $(TMPDIR)/WCSimRecoNtuple.o $(TMPDIR)/WCSimVertexNtuple.o $(TMPDIR)/WCSimVertexSeedNtuple.o $(TMPDIR)/WCSimNtupleWriter.o $(TMPDIR)/WCSimMsg.o  $(TMPDIR)/WCSimAnalysisRootDict.o 

ROOTEXTOBJS := $(WCSIM_LIBDIR)/WCSimRootEvent.o $(WCSIM_LIBDIR)/WCSimRootGeom.o $(WCSIM_LIBDIR)/WCSimRootDict.o

$(TMPDIR)/%.o : $(SRCDIR)/%.cc
	@echo "<**Compiling $@**>"
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) -c $< -o $@

$(TMPDIR)/%.d: $(SRCDIR)/%.cc
	@echo "<**Depend $@**>"
	set -e; $(CXX) $(CXXDEPEND) $(CXXFLAGS) $(INCLUDES) $(WCSIM_INCLUDES) $< \
	| sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

$(ROOTDICT) : $(ROOTSRC)
	rootcint -f $(ROOTDICT) -c -I$(INCDIR) -I$(WCSIM_INCDIR) -I$(shell root-config --incdir) WCSimIntegralLookupMaker3D.hh WCSimIntegralLookup3D.hh WCSimIntegralLookupReader.hh WCSimIntegralLookup.hh WCSimIntegralLookupMaker.hh WCSimRecoEvDisplay.hh WCSimRecoSummary.hh WCSimEmissionProfiles.hh WCSimFitterConfig.hh WCSimFitterInterface.hh WCSimFitterParameters.hh WCSimFitterTree.hh WCSimFitterPlots.hh WCSimTimeLikelihood.hh WCSimAnalysisConfig.hh WCSimDigitizerLikelihood.hh WCSimTotalLikelihood.hh WCSimLikelihoodTrack.hh WCSimLikelihoodDigit.hh WCSimLikelihoodDigitArray.hh WCSimChargePredictor.hh WCSimLikelihoodTuner.hh WCSimLikelihoodFitter.hh WCSimDisplayViewer.hh WCSimDisplayFactory.hh WCSimDisplay.hh WCSimDisplayAB.hh WCSimEveDisplay.hh WCSimEventWriter.hh WCSimGeometry.hh WCSimInterface.hh WCSimParameters.hh WCSimRecoObjectTable.hh WCSimRecoFactory.hh WCSimReco.hh WCSimRecoAB.hh WCSimRecoDigit.hh WCSimRecoCluster.hh WCSimRecoClusterDigit.hh WCSimRecoRing.hh WCSimRecoVertex.hh WCSimRecoEvent.hh WCSimTrueEvent.hh WCSimTrueTrack.hh WCSimHoughTransform.hh WCSimHoughTransformArray.hh WCSimDataCleaner.hh WCSimVertexFinder.hh WCSimVertexGeometry.hh WCSimVertexViewer.hh WCSimRingFinder.hh WCSimRingViewer.hh WCSimNtupleFactory.hh WCSimNtuple.hh WCSimRecoNtuple.hh WCSimVertexNtuple.hh WCSimVertexSeedNtuple.hh WCSimNtupleWriter.hh WCSimMsg.hh WCSimAnalysisRootLinkDef.hh

rootcint: $(ROOTDICT)

shared: $(ROOTDICT) $(ROOTSRC) $(ROOTOBJS)
	@echo "<**Shared**>"
	g++ -shared $(ROOTLIBS) $(ROOTGLIBS) -O $(ROOTOBJS) -o $(ROOTSO)

libWCSimAnalysis.a : $(ROOTOBJS)
	$(RM) $@
	ar clq $@ $(ROOTOBJS)
	
evDisp:	
	g++ `root-config --cflags --glibs --libs --evelibs` -I./include ${WCSIM_INCLUDES} -L${WCSIMANAHOME} -L${WCSIMHOME} -o evDisplay evDisplay.cc ${WCSIMHOME}/src/WCSimRootDict.cc ${WCSIMANAHOME}/src/WCSimAnalysisRootDict.cc -lWCSim -lWCSimAnalysis -lEG -lSpectrum -lMinuit

fitterProfile:
	g++ `root-config --cflags --glibs --libs --evelibs` -I./include ${WCSIM_INCLUDES} -L${WCSIMANAHOME} -L${WCSIMHOME} -L/home/ajperch/software/gperftools/lib/ -o fitterProfile runFitterProfile.cc ${WCSIMHOME}/src/WCSimRootDict.cc ${WCSIMANAHOME}/src/WCSimAnalysisRootDict.cc -lWCSim -lWCSimAnalysis -lEG -lSpectrum -lMinuit -lprofiler

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
