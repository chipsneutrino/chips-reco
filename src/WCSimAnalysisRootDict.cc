// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdIsrcdIWCSimAnalysisRootDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TClassAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "WCSimDigitizerLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimLikelihoodTrack.hh"
#include "WCSimLikelihoodDigit.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimChargeLikelihood.hh"
#include "WCSimLikelihoodTuner.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimDisplayViewer.hh"
#include "WCSimDisplayFactory.hh"
#include "WCSimDisplay.hh"
#include "WCSimDisplayAB.hh"
#include "WCSimEveDisplay.hh"
#include "WCSimEventWriter.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimParameters.hh"
#include "WCSimRecoObjectTable.hh"
#include "WCSimRecoFactory.hh"
#include "WCSimReco.hh"
#include "WCSimRecoAB.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoCluster.hh"
#include "WCSimRecoClusterDigit.hh"
#include "WCSimRecoRing.hh"
#include "WCSimRecoVertex.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimTrueEvent.hh"
#include "WCSimTrueTrack.hh"
#include "WCSimHoughTransform.hh"
#include "WCSimHoughTransformArray.hh"
#include "WCSimDataCleaner.hh"
#include "WCSimVertexFinder.hh"
#include "WCSimVertexGeometry.hh"
#include "WCSimVertexViewer.hh"
#include "WCSimRingFinder.hh"
#include "WCSimRingViewer.hh"
#include "WCSimNtupleFactory.hh"
#include "WCSimNtuple.hh"
#include "WCSimRecoNtuple.hh"
#include "WCSimVertexNtuple.hh"
#include "WCSimVertexSeedNtuple.hh"
#include "WCSimNtupleWriter.hh"
#include "WCSimMsg.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   void WCSimDisplay_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimDisplay(void *p = 0);
   static void *newArray_WCSimDisplay(Long_t size, void *p);
   static void delete_WCSimDisplay(void *p);
   static void deleteArray_WCSimDisplay(void *p);
   static void destruct_WCSimDisplay(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimDisplay*)
   {
      ::WCSimDisplay *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimDisplay >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimDisplay", ::WCSimDisplay::Class_Version(), "WCSimDisplay.hh", 9,
                  typeid(::WCSimDisplay), DefineBehavior(ptr, ptr),
                  &::WCSimDisplay::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimDisplay) );
      instance.SetNew(&new_WCSimDisplay);
      instance.SetNewArray(&newArray_WCSimDisplay);
      instance.SetDelete(&delete_WCSimDisplay);
      instance.SetDeleteArray(&deleteArray_WCSimDisplay);
      instance.SetDestructor(&destruct_WCSimDisplay);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimDisplay*)
   {
      return GenerateInitInstanceLocal((::WCSimDisplay*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimDisplay*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimDisplayFactory_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimDisplayFactory*)
   {
      ::WCSimDisplayFactory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimDisplayFactory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimDisplayFactory", ::WCSimDisplayFactory::Class_Version(), "WCSimDisplayFactory.hh", 8,
                  typeid(::WCSimDisplayFactory), DefineBehavior(ptr, ptr),
                  &::WCSimDisplayFactory::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimDisplayFactory) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimDisplayFactory*)
   {
      return GenerateInitInstanceLocal((::WCSimDisplayFactory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimDisplayFactory*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimDisplayViewer_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimDisplayViewer(void *p = 0);
   static void *newArray_WCSimDisplayViewer(Long_t size, void *p);
   static void delete_WCSimDisplayViewer(void *p);
   static void deleteArray_WCSimDisplayViewer(void *p);
   static void destruct_WCSimDisplayViewer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimDisplayViewer*)
   {
      ::WCSimDisplayViewer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimDisplayViewer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimDisplayViewer", ::WCSimDisplayViewer::Class_Version(), "WCSimDisplayViewer.hh", 24,
                  typeid(::WCSimDisplayViewer), DefineBehavior(ptr, ptr),
                  &::WCSimDisplayViewer::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimDisplayViewer) );
      instance.SetNew(&new_WCSimDisplayViewer);
      instance.SetNewArray(&newArray_WCSimDisplayViewer);
      instance.SetDelete(&delete_WCSimDisplayViewer);
      instance.SetDeleteArray(&deleteArray_WCSimDisplayViewer);
      instance.SetDestructor(&destruct_WCSimDisplayViewer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimDisplayViewer*)
   {
      return GenerateInitInstanceLocal((::WCSimDisplayViewer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimDisplayViewer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimDisplayAB_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimDisplayAB(void *p = 0);
   static void *newArray_WCSimDisplayAB(Long_t size, void *p);
   static void delete_WCSimDisplayAB(void *p);
   static void deleteArray_WCSimDisplayAB(void *p);
   static void destruct_WCSimDisplayAB(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimDisplayAB*)
   {
      ::WCSimDisplayAB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimDisplayAB >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimDisplayAB", ::WCSimDisplayAB::Class_Version(), "WCSimDisplayAB.hh", 24,
                  typeid(::WCSimDisplayAB), DefineBehavior(ptr, ptr),
                  &::WCSimDisplayAB::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimDisplayAB) );
      instance.SetNew(&new_WCSimDisplayAB);
      instance.SetNewArray(&newArray_WCSimDisplayAB);
      instance.SetDelete(&delete_WCSimDisplayAB);
      instance.SetDeleteArray(&deleteArray_WCSimDisplayAB);
      instance.SetDestructor(&destruct_WCSimDisplayAB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimDisplayAB*)
   {
      return GenerateInitInstanceLocal((::WCSimDisplayAB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimDisplayAB*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimEveDisplay_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimEveDisplay(void *p = 0);
   static void *newArray_WCSimEveDisplay(Long_t size, void *p);
   static void delete_WCSimEveDisplay(void *p);
   static void deleteArray_WCSimEveDisplay(void *p);
   static void destruct_WCSimEveDisplay(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimEveDisplay*)
   {
      ::WCSimEveDisplay *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimEveDisplay >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimEveDisplay", ::WCSimEveDisplay::Class_Version(), "WCSimEveDisplay.hh", 6,
                  typeid(::WCSimEveDisplay), DefineBehavior(ptr, ptr),
                  &::WCSimEveDisplay::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimEveDisplay) );
      instance.SetNew(&new_WCSimEveDisplay);
      instance.SetNewArray(&newArray_WCSimEveDisplay);
      instance.SetDelete(&delete_WCSimEveDisplay);
      instance.SetDeleteArray(&deleteArray_WCSimEveDisplay);
      instance.SetDestructor(&destruct_WCSimEveDisplay);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimEveDisplay*)
   {
      return GenerateInitInstanceLocal((::WCSimEveDisplay*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimEveDisplay*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimVertexViewer_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimVertexViewer(void *p = 0);
   static void *newArray_WCSimVertexViewer(Long_t size, void *p);
   static void delete_WCSimVertexViewer(void *p);
   static void deleteArray_WCSimVertexViewer(void *p);
   static void destruct_WCSimVertexViewer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimVertexViewer*)
   {
      ::WCSimVertexViewer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimVertexViewer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimVertexViewer", ::WCSimVertexViewer::Class_Version(), "WCSimVertexViewer.hh", 23,
                  typeid(::WCSimVertexViewer), DefineBehavior(ptr, ptr),
                  &::WCSimVertexViewer::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimVertexViewer) );
      instance.SetNew(&new_WCSimVertexViewer);
      instance.SetNewArray(&newArray_WCSimVertexViewer);
      instance.SetDelete(&delete_WCSimVertexViewer);
      instance.SetDeleteArray(&deleteArray_WCSimVertexViewer);
      instance.SetDestructor(&destruct_WCSimVertexViewer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimVertexViewer*)
   {
      return GenerateInitInstanceLocal((::WCSimVertexViewer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimVertexViewer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRingViewer_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimRingViewer(void *p = 0);
   static void *newArray_WCSimRingViewer(Long_t size, void *p);
   static void delete_WCSimRingViewer(void *p);
   static void deleteArray_WCSimRingViewer(void *p);
   static void destruct_WCSimRingViewer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRingViewer*)
   {
      ::WCSimRingViewer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRingViewer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRingViewer", ::WCSimRingViewer::Class_Version(), "WCSimRingViewer.hh", 15,
                  typeid(::WCSimRingViewer), DefineBehavior(ptr, ptr),
                  &::WCSimRingViewer::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRingViewer) );
      instance.SetNew(&new_WCSimRingViewer);
      instance.SetNewArray(&newArray_WCSimRingViewer);
      instance.SetDelete(&delete_WCSimRingViewer);
      instance.SetDeleteArray(&deleteArray_WCSimRingViewer);
      instance.SetDestructor(&destruct_WCSimRingViewer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRingViewer*)
   {
      return GenerateInitInstanceLocal((::WCSimRingViewer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRingViewer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimEventWriter_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimEventWriter(void *p = 0);
   static void *newArray_WCSimEventWriter(Long_t size, void *p);
   static void delete_WCSimEventWriter(void *p);
   static void deleteArray_WCSimEventWriter(void *p);
   static void destruct_WCSimEventWriter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimEventWriter*)
   {
      ::WCSimEventWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimEventWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimEventWriter", ::WCSimEventWriter::Class_Version(), "WCSimEventWriter.hh", 13,
                  typeid(::WCSimEventWriter), DefineBehavior(ptr, ptr),
                  &::WCSimEventWriter::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimEventWriter) );
      instance.SetNew(&new_WCSimEventWriter);
      instance.SetNewArray(&newArray_WCSimEventWriter);
      instance.SetDelete(&delete_WCSimEventWriter);
      instance.SetDeleteArray(&deleteArray_WCSimEventWriter);
      instance.SetDestructor(&destruct_WCSimEventWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimEventWriter*)
   {
      return GenerateInitInstanceLocal((::WCSimEventWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimEventWriter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimGeometry_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimGeometry*)
   {
      ::WCSimGeometry *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimGeometry >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimGeometry", ::WCSimGeometry::Class_Version(), "WCSimGeometry.hh", 8,
                  typeid(::WCSimGeometry), DefineBehavior(ptr, ptr),
                  &::WCSimGeometry::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimGeometry) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimGeometry*)
   {
      return GenerateInitInstanceLocal((::WCSimGeometry*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimGeometry*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimInterface_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimInterface*)
   {
      ::WCSimInterface *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimInterface >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimInterface", ::WCSimInterface::Class_Version(), "WCSimInterface.hh", 18,
                  typeid(::WCSimInterface), DefineBehavior(ptr, ptr),
                  &::WCSimInterface::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimInterface) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimInterface*)
   {
      return GenerateInitInstanceLocal((::WCSimInterface*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimInterface*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimParameters_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimParameters*)
   {
      ::WCSimParameters *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimParameters >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimParameters", ::WCSimParameters::Class_Version(), "WCSimParameters.hh", 6,
                  typeid(::WCSimParameters), DefineBehavior(ptr, ptr),
                  &::WCSimParameters::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimParameters) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimParameters*)
   {
      return GenerateInitInstanceLocal((::WCSimParameters*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimParameters*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimTrueEvent_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimTrueEvent(void *p = 0);
   static void *newArray_WCSimTrueEvent(Long_t size, void *p);
   static void delete_WCSimTrueEvent(void *p);
   static void deleteArray_WCSimTrueEvent(void *p);
   static void destruct_WCSimTrueEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTrueEvent*)
   {
      ::WCSimTrueEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimTrueEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTrueEvent", ::WCSimTrueEvent::Class_Version(), "WCSimTrueEvent.hh", 10,
                  typeid(::WCSimTrueEvent), DefineBehavior(ptr, ptr),
                  &::WCSimTrueEvent::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTrueEvent) );
      instance.SetNew(&new_WCSimTrueEvent);
      instance.SetNewArray(&newArray_WCSimTrueEvent);
      instance.SetDelete(&delete_WCSimTrueEvent);
      instance.SetDeleteArray(&deleteArray_WCSimTrueEvent);
      instance.SetDestructor(&destruct_WCSimTrueEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTrueEvent*)
   {
      return GenerateInitInstanceLocal((::WCSimTrueEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTrueEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimTrueTrack_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimTrueTrack(void *p);
   static void deleteArray_WCSimTrueTrack(void *p);
   static void destruct_WCSimTrueTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTrueTrack*)
   {
      ::WCSimTrueTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimTrueTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTrueTrack", ::WCSimTrueTrack::Class_Version(), "WCSimTrueTrack.hh", 6,
                  typeid(::WCSimTrueTrack), DefineBehavior(ptr, ptr),
                  &::WCSimTrueTrack::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTrueTrack) );
      instance.SetDelete(&delete_WCSimTrueTrack);
      instance.SetDeleteArray(&deleteArray_WCSimTrueTrack);
      instance.SetDestructor(&destruct_WCSimTrueTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTrueTrack*)
   {
      return GenerateInitInstanceLocal((::WCSimTrueTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTrueTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoObjectTable_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoObjectTable*)
   {
      ::WCSimRecoObjectTable *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoObjectTable >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoObjectTable", ::WCSimRecoObjectTable::Class_Version(), "WCSimRecoObjectTable.hh", 6,
                  typeid(::WCSimRecoObjectTable), DefineBehavior(ptr, ptr),
                  &::WCSimRecoObjectTable::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoObjectTable) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoObjectTable*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoObjectTable*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoFactory_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoFactory*)
   {
      ::WCSimRecoFactory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoFactory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoFactory", ::WCSimRecoFactory::Class_Version(), "WCSimRecoFactory.hh", 8,
                  typeid(::WCSimRecoFactory), DefineBehavior(ptr, ptr),
                  &::WCSimRecoFactory::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoFactory) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoFactory*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoFactory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoFactory*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimReco_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimReco(void *p);
   static void deleteArray_WCSimReco(void *p);
   static void destruct_WCSimReco(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimReco*)
   {
      ::WCSimReco *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimReco >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimReco", ::WCSimReco::Class_Version(), "WCSimReco.hh", 11,
                  typeid(::WCSimReco), DefineBehavior(ptr, ptr),
                  &::WCSimReco::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimReco) );
      instance.SetDelete(&delete_WCSimReco);
      instance.SetDeleteArray(&deleteArray_WCSimReco);
      instance.SetDestructor(&destruct_WCSimReco);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimReco*)
   {
      return GenerateInitInstanceLocal((::WCSimReco*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimReco*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoAB_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimRecoAB(void *p = 0);
   static void *newArray_WCSimRecoAB(Long_t size, void *p);
   static void delete_WCSimRecoAB(void *p);
   static void deleteArray_WCSimRecoAB(void *p);
   static void destruct_WCSimRecoAB(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoAB*)
   {
      ::WCSimRecoAB *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoAB >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoAB", ::WCSimRecoAB::Class_Version(), "WCSimRecoAB.hh", 12,
                  typeid(::WCSimRecoAB), DefineBehavior(ptr, ptr),
                  &::WCSimRecoAB::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoAB) );
      instance.SetNew(&new_WCSimRecoAB);
      instance.SetNewArray(&newArray_WCSimRecoAB);
      instance.SetDelete(&delete_WCSimRecoAB);
      instance.SetDeleteArray(&deleteArray_WCSimRecoAB);
      instance.SetDestructor(&destruct_WCSimRecoAB);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoAB*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoAB*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoAB*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoEvent_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimRecoEvent(void *p = 0);
   static void *newArray_WCSimRecoEvent(Long_t size, void *p);
   static void delete_WCSimRecoEvent(void *p);
   static void deleteArray_WCSimRecoEvent(void *p);
   static void destruct_WCSimRecoEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoEvent*)
   {
      ::WCSimRecoEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoEvent", ::WCSimRecoEvent::Class_Version(), "WCSimRecoEvent.hh", 12,
                  typeid(::WCSimRecoEvent), DefineBehavior(ptr, ptr),
                  &::WCSimRecoEvent::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoEvent) );
      instance.SetNew(&new_WCSimRecoEvent);
      instance.SetNewArray(&newArray_WCSimRecoEvent);
      instance.SetDelete(&delete_WCSimRecoEvent);
      instance.SetDeleteArray(&deleteArray_WCSimRecoEvent);
      instance.SetDestructor(&destruct_WCSimRecoEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoEvent*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoDigit_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimRecoDigit(void *p);
   static void deleteArray_WCSimRecoDigit(void *p);
   static void destruct_WCSimRecoDigit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoDigit*)
   {
      ::WCSimRecoDigit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoDigit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoDigit", ::WCSimRecoDigit::Class_Version(), "WCSimRecoDigit.hh", 6,
                  typeid(::WCSimRecoDigit), DefineBehavior(ptr, ptr),
                  &::WCSimRecoDigit::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoDigit) );
      instance.SetDelete(&delete_WCSimRecoDigit);
      instance.SetDeleteArray(&deleteArray_WCSimRecoDigit);
      instance.SetDestructor(&destruct_WCSimRecoDigit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoDigit*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoDigit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoCluster_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimRecoCluster(void *p = 0);
   static void *newArray_WCSimRecoCluster(Long_t size, void *p);
   static void delete_WCSimRecoCluster(void *p);
   static void deleteArray_WCSimRecoCluster(void *p);
   static void destruct_WCSimRecoCluster(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoCluster*)
   {
      ::WCSimRecoCluster *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoCluster >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoCluster", ::WCSimRecoCluster::Class_Version(), "WCSimRecoCluster.hh", 10,
                  typeid(::WCSimRecoCluster), DefineBehavior(ptr, ptr),
                  &::WCSimRecoCluster::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoCluster) );
      instance.SetNew(&new_WCSimRecoCluster);
      instance.SetNewArray(&newArray_WCSimRecoCluster);
      instance.SetDelete(&delete_WCSimRecoCluster);
      instance.SetDeleteArray(&deleteArray_WCSimRecoCluster);
      instance.SetDestructor(&destruct_WCSimRecoCluster);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoCluster*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoCluster*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoCluster*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoClusterDigit_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimRecoClusterDigit(void *p);
   static void deleteArray_WCSimRecoClusterDigit(void *p);
   static void destruct_WCSimRecoClusterDigit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoClusterDigit*)
   {
      ::WCSimRecoClusterDigit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoClusterDigit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoClusterDigit", ::WCSimRecoClusterDigit::Class_Version(), "WCSimRecoClusterDigit.hh", 10,
                  typeid(::WCSimRecoClusterDigit), DefineBehavior(ptr, ptr),
                  &::WCSimRecoClusterDigit::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoClusterDigit) );
      instance.SetDelete(&delete_WCSimRecoClusterDigit);
      instance.SetDeleteArray(&deleteArray_WCSimRecoClusterDigit);
      instance.SetDestructor(&destruct_WCSimRecoClusterDigit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoClusterDigit*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoClusterDigit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoClusterDigit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoRing_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimRecoRing(void *p);
   static void deleteArray_WCSimRecoRing(void *p);
   static void destruct_WCSimRecoRing(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoRing*)
   {
      ::WCSimRecoRing *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoRing >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoRing", ::WCSimRecoRing::Class_Version(), "WCSimRecoRing.hh", 6,
                  typeid(::WCSimRecoRing), DefineBehavior(ptr, ptr),
                  &::WCSimRecoRing::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoRing) );
      instance.SetDelete(&delete_WCSimRecoRing);
      instance.SetDeleteArray(&deleteArray_WCSimRecoRing);
      instance.SetDestructor(&destruct_WCSimRecoRing);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoRing*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoRing*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoRing*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoVertex_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimRecoVertex(void *p = 0);
   static void *newArray_WCSimRecoVertex(Long_t size, void *p);
   static void delete_WCSimRecoVertex(void *p);
   static void deleteArray_WCSimRecoVertex(void *p);
   static void destruct_WCSimRecoVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoVertex*)
   {
      ::WCSimRecoVertex *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoVertex >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoVertex", ::WCSimRecoVertex::Class_Version(), "WCSimRecoVertex.hh", 6,
                  typeid(::WCSimRecoVertex), DefineBehavior(ptr, ptr),
                  &::WCSimRecoVertex::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoVertex) );
      instance.SetNew(&new_WCSimRecoVertex);
      instance.SetNewArray(&newArray_WCSimRecoVertex);
      instance.SetDelete(&delete_WCSimRecoVertex);
      instance.SetDeleteArray(&deleteArray_WCSimRecoVertex);
      instance.SetDestructor(&destruct_WCSimRecoVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoVertex*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoVertex*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoVertex*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimVertexFinder_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimVertexFinder*)
   {
      ::WCSimVertexFinder *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimVertexFinder >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimVertexFinder", ::WCSimVertexFinder::Class_Version(), "WCSimVertexFinder.hh", 14,
                  typeid(::WCSimVertexFinder), DefineBehavior(ptr, ptr),
                  &::WCSimVertexFinder::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimVertexFinder) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimVertexFinder*)
   {
      return GenerateInitInstanceLocal((::WCSimVertexFinder*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimVertexFinder*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRingFinder_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRingFinder*)
   {
      ::WCSimRingFinder *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRingFinder >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRingFinder", ::WCSimRingFinder::Class_Version(), "WCSimRingFinder.hh", 15,
                  typeid(::WCSimRingFinder), DefineBehavior(ptr, ptr),
                  &::WCSimRingFinder::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRingFinder) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRingFinder*)
   {
      return GenerateInitInstanceLocal((::WCSimRingFinder*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRingFinder*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimDataCleaner_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimDataCleaner*)
   {
      ::WCSimDataCleaner *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimDataCleaner >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimDataCleaner", ::WCSimDataCleaner::Class_Version(), "WCSimDataCleaner.hh", 12,
                  typeid(::WCSimDataCleaner), DefineBehavior(ptr, ptr),
                  &::WCSimDataCleaner::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimDataCleaner) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimDataCleaner*)
   {
      return GenerateInitInstanceLocal((::WCSimDataCleaner*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimDataCleaner*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimVertexGeometry_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimVertexGeometry*)
   {
      ::WCSimVertexGeometry *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimVertexGeometry >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimVertexGeometry", ::WCSimVertexGeometry::Class_Version(), "WCSimVertexGeometry.hh", 11,
                  typeid(::WCSimVertexGeometry), DefineBehavior(ptr, ptr),
                  &::WCSimVertexGeometry::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimVertexGeometry) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimVertexGeometry*)
   {
      return GenerateInitInstanceLocal((::WCSimVertexGeometry*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimVertexGeometry*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimHoughTransform_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimHoughTransform(void *p);
   static void deleteArray_WCSimHoughTransform(void *p);
   static void destruct_WCSimHoughTransform(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimHoughTransform*)
   {
      ::WCSimHoughTransform *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimHoughTransform >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimHoughTransform", ::WCSimHoughTransform::Class_Version(), "WCSimHoughTransform.hh", 8,
                  typeid(::WCSimHoughTransform), DefineBehavior(ptr, ptr),
                  &::WCSimHoughTransform::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimHoughTransform) );
      instance.SetDelete(&delete_WCSimHoughTransform);
      instance.SetDeleteArray(&deleteArray_WCSimHoughTransform);
      instance.SetDestructor(&destruct_WCSimHoughTransform);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimHoughTransform*)
   {
      return GenerateInitInstanceLocal((::WCSimHoughTransform*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimHoughTransform*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimHoughTransformArray_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimHoughTransformArray(void *p);
   static void deleteArray_WCSimHoughTransformArray(void *p);
   static void destruct_WCSimHoughTransformArray(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimHoughTransformArray*)
   {
      ::WCSimHoughTransformArray *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimHoughTransformArray >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimHoughTransformArray", ::WCSimHoughTransformArray::Class_Version(), "WCSimHoughTransformArray.hh", 10,
                  typeid(::WCSimHoughTransformArray), DefineBehavior(ptr, ptr),
                  &::WCSimHoughTransformArray::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimHoughTransformArray) );
      instance.SetDelete(&delete_WCSimHoughTransformArray);
      instance.SetDeleteArray(&deleteArray_WCSimHoughTransformArray);
      instance.SetDestructor(&destruct_WCSimHoughTransformArray);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimHoughTransformArray*)
   {
      return GenerateInitInstanceLocal((::WCSimHoughTransformArray*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimHoughTransformArray*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimNtupleFactory_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimNtupleFactory*)
   {
      ::WCSimNtupleFactory *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimNtupleFactory >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimNtupleFactory", ::WCSimNtupleFactory::Class_Version(), "WCSimNtupleFactory.hh", 8,
                  typeid(::WCSimNtupleFactory), DefineBehavior(ptr, ptr),
                  &::WCSimNtupleFactory::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimNtupleFactory) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimNtupleFactory*)
   {
      return GenerateInitInstanceLocal((::WCSimNtupleFactory*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimNtupleFactory*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimNtupleWriter_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimNtupleWriter(void *p = 0);
   static void *newArray_WCSimNtupleWriter(Long_t size, void *p);
   static void delete_WCSimNtupleWriter(void *p);
   static void deleteArray_WCSimNtupleWriter(void *p);
   static void destruct_WCSimNtupleWriter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimNtupleWriter*)
   {
      ::WCSimNtupleWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimNtupleWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimNtupleWriter", ::WCSimNtupleWriter::Class_Version(), "WCSimNtupleWriter.hh", 13,
                  typeid(::WCSimNtupleWriter), DefineBehavior(ptr, ptr),
                  &::WCSimNtupleWriter::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimNtupleWriter) );
      instance.SetNew(&new_WCSimNtupleWriter);
      instance.SetNewArray(&newArray_WCSimNtupleWriter);
      instance.SetDelete(&delete_WCSimNtupleWriter);
      instance.SetDeleteArray(&deleteArray_WCSimNtupleWriter);
      instance.SetDestructor(&destruct_WCSimNtupleWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimNtupleWriter*)
   {
      return GenerateInitInstanceLocal((::WCSimNtupleWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimNtupleWriter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimNtuple_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimNtuple(void *p = 0);
   static void *newArray_WCSimNtuple(Long_t size, void *p);
   static void delete_WCSimNtuple(void *p);
   static void deleteArray_WCSimNtuple(void *p);
   static void destruct_WCSimNtuple(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimNtuple*)
   {
      ::WCSimNtuple *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimNtuple >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimNtuple", ::WCSimNtuple::Class_Version(), "WCSimNtuple.hh", 10,
                  typeid(::WCSimNtuple), DefineBehavior(ptr, ptr),
                  &::WCSimNtuple::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimNtuple) );
      instance.SetNew(&new_WCSimNtuple);
      instance.SetNewArray(&newArray_WCSimNtuple);
      instance.SetDelete(&delete_WCSimNtuple);
      instance.SetDeleteArray(&deleteArray_WCSimNtuple);
      instance.SetDestructor(&destruct_WCSimNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimNtuple*)
   {
      return GenerateInitInstanceLocal((::WCSimNtuple*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimNtuple*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimRecoNtuple_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimRecoNtuple(void *p = 0);
   static void *newArray_WCSimRecoNtuple(Long_t size, void *p);
   static void delete_WCSimRecoNtuple(void *p);
   static void deleteArray_WCSimRecoNtuple(void *p);
   static void destruct_WCSimRecoNtuple(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoNtuple*)
   {
      ::WCSimRecoNtuple *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoNtuple >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoNtuple", ::WCSimRecoNtuple::Class_Version(), "WCSimRecoNtuple.hh", 12,
                  typeid(::WCSimRecoNtuple), DefineBehavior(ptr, ptr),
                  &::WCSimRecoNtuple::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoNtuple) );
      instance.SetNew(&new_WCSimRecoNtuple);
      instance.SetNewArray(&newArray_WCSimRecoNtuple);
      instance.SetDelete(&delete_WCSimRecoNtuple);
      instance.SetDeleteArray(&deleteArray_WCSimRecoNtuple);
      instance.SetDestructor(&destruct_WCSimRecoNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoNtuple*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoNtuple*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoNtuple*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimVertexNtuple_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimVertexNtuple(void *p = 0);
   static void *newArray_WCSimVertexNtuple(Long_t size, void *p);
   static void delete_WCSimVertexNtuple(void *p);
   static void deleteArray_WCSimVertexNtuple(void *p);
   static void destruct_WCSimVertexNtuple(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimVertexNtuple*)
   {
      ::WCSimVertexNtuple *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimVertexNtuple >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimVertexNtuple", ::WCSimVertexNtuple::Class_Version(), "WCSimVertexNtuple.hh", 12,
                  typeid(::WCSimVertexNtuple), DefineBehavior(ptr, ptr),
                  &::WCSimVertexNtuple::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimVertexNtuple) );
      instance.SetNew(&new_WCSimVertexNtuple);
      instance.SetNewArray(&newArray_WCSimVertexNtuple);
      instance.SetDelete(&delete_WCSimVertexNtuple);
      instance.SetDeleteArray(&deleteArray_WCSimVertexNtuple);
      instance.SetDestructor(&destruct_WCSimVertexNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimVertexNtuple*)
   {
      return GenerateInitInstanceLocal((::WCSimVertexNtuple*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimVertexNtuple*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimVertexSeedNtuple_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimVertexSeedNtuple(void *p = 0);
   static void *newArray_WCSimVertexSeedNtuple(Long_t size, void *p);
   static void delete_WCSimVertexSeedNtuple(void *p);
   static void deleteArray_WCSimVertexSeedNtuple(void *p);
   static void destruct_WCSimVertexSeedNtuple(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimVertexSeedNtuple*)
   {
      ::WCSimVertexSeedNtuple *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimVertexSeedNtuple >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimVertexSeedNtuple", ::WCSimVertexSeedNtuple::Class_Version(), "WCSimVertexSeedNtuple.hh", 12,
                  typeid(::WCSimVertexSeedNtuple), DefineBehavior(ptr, ptr),
                  &::WCSimVertexSeedNtuple::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimVertexSeedNtuple) );
      instance.SetNew(&new_WCSimVertexSeedNtuple);
      instance.SetNewArray(&newArray_WCSimVertexSeedNtuple);
      instance.SetDelete(&delete_WCSimVertexSeedNtuple);
      instance.SetDeleteArray(&deleteArray_WCSimVertexSeedNtuple);
      instance.SetDestructor(&destruct_WCSimVertexSeedNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimVertexSeedNtuple*)
   {
      return GenerateInitInstanceLocal((::WCSimVertexSeedNtuple*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimVertexSeedNtuple*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimMsg_ShowMembers(void *obj, TMemberInspector &R__insp);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimMsg*)
   {
      ::WCSimMsg *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimMsg >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimMsg", ::WCSimMsg::Class_Version(), "WCSimMsg.hh", 6,
                  typeid(::WCSimMsg), DefineBehavior(ptr, ptr),
                  &::WCSimMsg::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimMsg) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimMsg*)
   {
      return GenerateInitInstanceLocal((::WCSimMsg*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimMsg*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimTotalLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimTotalLikelihood(void *p);
   static void deleteArray_WCSimTotalLikelihood(void *p);
   static void destruct_WCSimTotalLikelihood(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTotalLikelihood*)
   {
      ::WCSimTotalLikelihood *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimTotalLikelihood >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTotalLikelihood", ::WCSimTotalLikelihood::Class_Version(), "WCSimTotalLikelihood.hh", 12,
                  typeid(::WCSimTotalLikelihood), DefineBehavior(ptr, ptr),
                  &::WCSimTotalLikelihood::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTotalLikelihood) );
      instance.SetDelete(&delete_WCSimTotalLikelihood);
      instance.SetDeleteArray(&deleteArray_WCSimTotalLikelihood);
      instance.SetDestructor(&destruct_WCSimTotalLikelihood);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTotalLikelihood*)
   {
      return GenerateInitInstanceLocal((::WCSimTotalLikelihood*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTotalLikelihood*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimLikelihoodFitter_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimLikelihoodFitter(void *p);
   static void deleteArray_WCSimLikelihoodFitter(void *p);
   static void destruct_WCSimLikelihoodFitter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimLikelihoodFitter*)
   {
      ::WCSimLikelihoodFitter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimLikelihoodFitter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimLikelihoodFitter", ::WCSimLikelihoodFitter::Class_Version(), "WCSimLikelihoodFitter.hh", 15,
                  typeid(::WCSimLikelihoodFitter), DefineBehavior(ptr, ptr),
                  &::WCSimLikelihoodFitter::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimLikelihoodFitter) );
      instance.SetDelete(&delete_WCSimLikelihoodFitter);
      instance.SetDeleteArray(&deleteArray_WCSimLikelihoodFitter);
      instance.SetDestructor(&destruct_WCSimLikelihoodFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimLikelihoodFitter*)
   {
      return GenerateInitInstanceLocal((::WCSimLikelihoodFitter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimLikelihoodFitter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimLikelihoodTuner_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void WCSimLikelihoodTuner_Dictionary();
   static void WCSimLikelihoodTuner_TClassManip(TClass*);
   static void delete_WCSimLikelihoodTuner(void *p);
   static void deleteArray_WCSimLikelihoodTuner(void *p);
   static void destruct_WCSimLikelihoodTuner(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimLikelihoodTuner*)
   {
      ::WCSimLikelihoodTuner *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WCSimLikelihoodTuner),0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimLikelihoodTuner", "WCSimLikelihoodTuner.hh", 17,
                  typeid(::WCSimLikelihoodTuner), DefineBehavior(ptr, ptr),
                  &WCSimLikelihoodTuner_ShowMembers, &WCSimLikelihoodTuner_Dictionary, isa_proxy, 4,
                  sizeof(::WCSimLikelihoodTuner) );
      instance.SetDelete(&delete_WCSimLikelihoodTuner);
      instance.SetDeleteArray(&deleteArray_WCSimLikelihoodTuner);
      instance.SetDestructor(&destruct_WCSimLikelihoodTuner);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimLikelihoodTuner*)
   {
      return GenerateInitInstanceLocal((::WCSimLikelihoodTuner*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimLikelihoodTuner*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void WCSimLikelihoodTuner_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodTuner*)0x0)->GetClass();
      WCSimLikelihoodTuner_TClassManip(theClass);
   }

   static void WCSimLikelihoodTuner_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   void WCSimChargeLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void WCSimChargeLikelihood_Dictionary();
   static void WCSimChargeLikelihood_TClassManip(TClass*);
   static void delete_WCSimChargeLikelihood(void *p);
   static void deleteArray_WCSimChargeLikelihood(void *p);
   static void destruct_WCSimChargeLikelihood(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimChargeLikelihood*)
   {
      ::WCSimChargeLikelihood *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WCSimChargeLikelihood),0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimChargeLikelihood", "WCSimChargeLikelihood.hh", 19,
                  typeid(::WCSimChargeLikelihood), DefineBehavior(ptr, ptr),
                  &WCSimChargeLikelihood_ShowMembers, &WCSimChargeLikelihood_Dictionary, isa_proxy, 4,
                  sizeof(::WCSimChargeLikelihood) );
      instance.SetDelete(&delete_WCSimChargeLikelihood);
      instance.SetDeleteArray(&deleteArray_WCSimChargeLikelihood);
      instance.SetDestructor(&destruct_WCSimChargeLikelihood);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimChargeLikelihood*)
   {
      return GenerateInitInstanceLocal((::WCSimChargeLikelihood*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimChargeLikelihood*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void WCSimChargeLikelihood_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WCSimChargeLikelihood*)0x0)->GetClass();
      WCSimChargeLikelihood_TClassManip(theClass);
   }

   static void WCSimChargeLikelihood_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   void WCSimLikelihoodDigit_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_WCSimLikelihoodDigit(void *p);
   static void deleteArray_WCSimLikelihoodDigit(void *p);
   static void destruct_WCSimLikelihoodDigit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimLikelihoodDigit*)
   {
      ::WCSimLikelihoodDigit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimLikelihoodDigit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimLikelihoodDigit", ::WCSimLikelihoodDigit::Class_Version(), "WCSimLikelihoodDigit.hh", 9,
                  typeid(::WCSimLikelihoodDigit), DefineBehavior(ptr, ptr),
                  &::WCSimLikelihoodDigit::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimLikelihoodDigit) );
      instance.SetDelete(&delete_WCSimLikelihoodDigit);
      instance.SetDeleteArray(&deleteArray_WCSimLikelihoodDigit);
      instance.SetDestructor(&destruct_WCSimLikelihoodDigit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimLikelihoodDigit*)
   {
      return GenerateInitInstanceLocal((::WCSimLikelihoodDigit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimLikelihoodDigit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimLikelihoodDigitArray_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimLikelihoodDigitArray(void *p = 0);
   static void *newArray_WCSimLikelihoodDigitArray(Long_t size, void *p);
   static void delete_WCSimLikelihoodDigitArray(void *p);
   static void deleteArray_WCSimLikelihoodDigitArray(void *p);
   static void destruct_WCSimLikelihoodDigitArray(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimLikelihoodDigitArray*)
   {
      ::WCSimLikelihoodDigitArray *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimLikelihoodDigitArray >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimLikelihoodDigitArray", ::WCSimLikelihoodDigitArray::Class_Version(), "WCSimLikelihoodDigitArray.hh", 10,
                  typeid(::WCSimLikelihoodDigitArray), DefineBehavior(ptr, ptr),
                  &::WCSimLikelihoodDigitArray::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimLikelihoodDigitArray) );
      instance.SetNew(&new_WCSimLikelihoodDigitArray);
      instance.SetNewArray(&newArray_WCSimLikelihoodDigitArray);
      instance.SetDelete(&delete_WCSimLikelihoodDigitArray);
      instance.SetDeleteArray(&deleteArray_WCSimLikelihoodDigitArray);
      instance.SetDestructor(&destruct_WCSimLikelihoodDigitArray);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimLikelihoodDigitArray*)
   {
      return GenerateInitInstanceLocal((::WCSimLikelihoodDigitArray*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimLikelihoodDigitArray*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimLikelihoodTrack_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_WCSimLikelihoodTrack(void *p = 0);
   static void *newArray_WCSimLikelihoodTrack(Long_t size, void *p);
   static void delete_WCSimLikelihoodTrack(void *p);
   static void deleteArray_WCSimLikelihoodTrack(void *p);
   static void destruct_WCSimLikelihoodTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimLikelihoodTrack*)
   {
      ::WCSimLikelihoodTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimLikelihoodTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimLikelihoodTrack", ::WCSimLikelihoodTrack::Class_Version(), "WCSimLikelihoodTrack.hh", 8,
                  typeid(::WCSimLikelihoodTrack), DefineBehavior(ptr, ptr),
                  &::WCSimLikelihoodTrack::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimLikelihoodTrack) );
      instance.SetNew(&new_WCSimLikelihoodTrack);
      instance.SetNewArray(&newArray_WCSimLikelihoodTrack);
      instance.SetDelete(&delete_WCSimLikelihoodTrack);
      instance.SetDeleteArray(&deleteArray_WCSimLikelihoodTrack);
      instance.SetDestructor(&destruct_WCSimLikelihoodTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimLikelihoodTrack*)
   {
      return GenerateInitInstanceLocal((::WCSimLikelihoodTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimLikelihoodTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   void WCSimTimeLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void WCSimTimeLikelihood_Dictionary();
   static void WCSimTimeLikelihood_TClassManip(TClass*);
   static void delete_WCSimTimeLikelihood(void *p);
   static void deleteArray_WCSimTimeLikelihood(void *p);
   static void destruct_WCSimTimeLikelihood(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTimeLikelihood*)
   {
      ::WCSimTimeLikelihood *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WCSimTimeLikelihood),0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTimeLikelihood", "WCSimTimeLikelihood.hh", 16,
                  typeid(::WCSimTimeLikelihood), DefineBehavior(ptr, ptr),
                  &WCSimTimeLikelihood_ShowMembers, &WCSimTimeLikelihood_Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTimeLikelihood) );
      instance.SetDelete(&delete_WCSimTimeLikelihood);
      instance.SetDeleteArray(&deleteArray_WCSimTimeLikelihood);
      instance.SetDestructor(&destruct_WCSimTimeLikelihood);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTimeLikelihood*)
   {
      return GenerateInitInstanceLocal((::WCSimTimeLikelihood*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTimeLikelihood*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void WCSimTimeLikelihood_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WCSimTimeLikelihood*)0x0)->GetClass();
      WCSimTimeLikelihood_TClassManip(theClass);
   }

   static void WCSimTimeLikelihood_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   void WCSimDigitizerLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void WCSimDigitizerLikelihood_Dictionary();
   static void WCSimDigitizerLikelihood_TClassManip(TClass*);
   static void *new_WCSimDigitizerLikelihood(void *p = 0);
   static void *newArray_WCSimDigitizerLikelihood(Long_t size, void *p);
   static void delete_WCSimDigitizerLikelihood(void *p);
   static void deleteArray_WCSimDigitizerLikelihood(void *p);
   static void destruct_WCSimDigitizerLikelihood(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimDigitizerLikelihood*)
   {
      ::WCSimDigitizerLikelihood *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::WCSimDigitizerLikelihood),0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimDigitizerLikelihood", "WCSimDigitizerLikelihood.hh", 8,
                  typeid(::WCSimDigitizerLikelihood), DefineBehavior(ptr, ptr),
                  &WCSimDigitizerLikelihood_ShowMembers, &WCSimDigitizerLikelihood_Dictionary, isa_proxy, 4,
                  sizeof(::WCSimDigitizerLikelihood) );
      instance.SetNew(&new_WCSimDigitizerLikelihood);
      instance.SetNewArray(&newArray_WCSimDigitizerLikelihood);
      instance.SetDelete(&delete_WCSimDigitizerLikelihood);
      instance.SetDeleteArray(&deleteArray_WCSimDigitizerLikelihood);
      instance.SetDestructor(&destruct_WCSimDigitizerLikelihood);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimDigitizerLikelihood*)
   {
      return GenerateInitInstanceLocal((::WCSimDigitizerLikelihood*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimDigitizerLikelihood*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void WCSimDigitizerLikelihood_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::WCSimDigitizerLikelihood*)0x0)->GetClass();
      WCSimDigitizerLikelihood_TClassManip(theClass);
   }

   static void WCSimDigitizerLikelihood_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
TClass *WCSimDisplay::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimDisplay::Class_Name()
{
   return "WCSimDisplay";
}

//______________________________________________________________________________
const char *WCSimDisplay::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplay*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimDisplay::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplay*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimDisplay::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplay*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimDisplay::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplay*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimDisplayFactory::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimDisplayFactory::Class_Name()
{
   return "WCSimDisplayFactory";
}

//______________________________________________________________________________
const char *WCSimDisplayFactory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayFactory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimDisplayFactory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayFactory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimDisplayFactory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayFactory*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimDisplayFactory::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayFactory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimDisplayViewer::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimDisplayViewer::Class_Name()
{
   return "WCSimDisplayViewer";
}

//______________________________________________________________________________
const char *WCSimDisplayViewer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayViewer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimDisplayViewer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayViewer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimDisplayViewer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayViewer*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimDisplayViewer::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayViewer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimDisplayAB::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimDisplayAB::Class_Name()
{
   return "WCSimDisplayAB";
}

//______________________________________________________________________________
const char *WCSimDisplayAB::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayAB*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimDisplayAB::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayAB*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimDisplayAB::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayAB*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimDisplayAB::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDisplayAB*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimEveDisplay::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimEveDisplay::Class_Name()
{
   return "WCSimEveDisplay";
}

//______________________________________________________________________________
const char *WCSimEveDisplay::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimEveDisplay*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimEveDisplay::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimEveDisplay*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimEveDisplay::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimEveDisplay*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimEveDisplay::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimEveDisplay*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimVertexViewer::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimVertexViewer::Class_Name()
{
   return "WCSimVertexViewer";
}

//______________________________________________________________________________
const char *WCSimVertexViewer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexViewer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimVertexViewer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexViewer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimVertexViewer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexViewer*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimVertexViewer::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexViewer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRingViewer::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRingViewer::Class_Name()
{
   return "WCSimRingViewer";
}

//______________________________________________________________________________
const char *WCSimRingViewer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingViewer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRingViewer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingViewer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRingViewer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingViewer*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRingViewer::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingViewer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimEventWriter::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimEventWriter::Class_Name()
{
   return "WCSimEventWriter";
}

//______________________________________________________________________________
const char *WCSimEventWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimEventWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimEventWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimEventWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimEventWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimEventWriter*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimEventWriter::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimEventWriter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimGeometry::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimGeometry::Class_Name()
{
   return "WCSimGeometry";
}

//______________________________________________________________________________
const char *WCSimGeometry::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimGeometry*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimGeometry::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimGeometry*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimGeometry::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimGeometry*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimGeometry::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimGeometry*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimInterface::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimInterface::Class_Name()
{
   return "WCSimInterface";
}

//______________________________________________________________________________
const char *WCSimInterface::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimInterface*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimInterface::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimInterface*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimInterface::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimInterface*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimInterface::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimInterface*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimParameters::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimParameters::Class_Name()
{
   return "WCSimParameters";
}

//______________________________________________________________________________
const char *WCSimParameters::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimParameters*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimParameters::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimParameters*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimParameters::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimParameters*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimParameters::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimParameters*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimTrueEvent::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimTrueEvent::Class_Name()
{
   return "WCSimTrueEvent";
}

//______________________________________________________________________________
const char *WCSimTrueEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimTrueEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimTrueEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueEvent*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimTrueEvent::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimTrueTrack::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimTrueTrack::Class_Name()
{
   return "WCSimTrueTrack";
}

//______________________________________________________________________________
const char *WCSimTrueTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimTrueTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimTrueTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueTrack*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimTrueTrack::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoObjectTable::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoObjectTable::Class_Name()
{
   return "WCSimRecoObjectTable";
}

//______________________________________________________________________________
const char *WCSimRecoObjectTable::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoObjectTable::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoObjectTable::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoObjectTable::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoFactory::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoFactory::Class_Name()
{
   return "WCSimRecoFactory";
}

//______________________________________________________________________________
const char *WCSimRecoFactory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoFactory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoFactory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoFactory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoFactory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoFactory*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoFactory::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoFactory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimReco::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimReco::Class_Name()
{
   return "WCSimReco";
}

//______________________________________________________________________________
const char *WCSimReco::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimReco*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimReco::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimReco*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimReco::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimReco*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimReco::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimReco*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoAB::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoAB::Class_Name()
{
   return "WCSimRecoAB";
}

//______________________________________________________________________________
const char *WCSimRecoAB::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoAB*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoAB::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoAB*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoAB::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoAB*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoAB::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoAB*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoEvent::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoEvent::Class_Name()
{
   return "WCSimRecoEvent";
}

//______________________________________________________________________________
const char *WCSimRecoEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoEvent*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoEvent::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoDigit::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoDigit::Class_Name()
{
   return "WCSimRecoDigit";
}

//______________________________________________________________________________
const char *WCSimRecoDigit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoDigit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoDigit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoDigit::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoCluster::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoCluster::Class_Name()
{
   return "WCSimRecoCluster";
}

//______________________________________________________________________________
const char *WCSimRecoCluster::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoCluster*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoCluster::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoCluster*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoCluster::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoCluster*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoCluster::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoCluster*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoClusterDigit::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoClusterDigit::Class_Name()
{
   return "WCSimRecoClusterDigit";
}

//______________________________________________________________________________
const char *WCSimRecoClusterDigit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoClusterDigit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoClusterDigit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoClusterDigit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoClusterDigit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoClusterDigit*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoClusterDigit::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoClusterDigit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoRing::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoRing::Class_Name()
{
   return "WCSimRecoRing";
}

//______________________________________________________________________________
const char *WCSimRecoRing::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoRing*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoRing::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoRing*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoRing::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoRing*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoRing::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoRing*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoVertex::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoVertex::Class_Name()
{
   return "WCSimRecoVertex";
}

//______________________________________________________________________________
const char *WCSimRecoVertex::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoVertex*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoVertex::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoVertex*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoVertex::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoVertex*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoVertex::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoVertex*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimVertexFinder::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimVertexFinder::Class_Name()
{
   return "WCSimVertexFinder";
}

//______________________________________________________________________________
const char *WCSimVertexFinder::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexFinder*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimVertexFinder::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexFinder*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimVertexFinder::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexFinder*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimVertexFinder::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexFinder*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRingFinder::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRingFinder::Class_Name()
{
   return "WCSimRingFinder";
}

//______________________________________________________________________________
const char *WCSimRingFinder::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingFinder*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRingFinder::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingFinder*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRingFinder::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingFinder*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRingFinder::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRingFinder*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimDataCleaner::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimDataCleaner::Class_Name()
{
   return "WCSimDataCleaner";
}

//______________________________________________________________________________
const char *WCSimDataCleaner::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDataCleaner*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimDataCleaner::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimDataCleaner*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimDataCleaner::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDataCleaner*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimDataCleaner::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimDataCleaner*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimVertexGeometry::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimVertexGeometry::Class_Name()
{
   return "WCSimVertexGeometry";
}

//______________________________________________________________________________
const char *WCSimVertexGeometry::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexGeometry*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimVertexGeometry::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexGeometry*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimVertexGeometry::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexGeometry*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimVertexGeometry::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexGeometry*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimHoughTransform::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimHoughTransform::Class_Name()
{
   return "WCSimHoughTransform";
}

//______________________________________________________________________________
const char *WCSimHoughTransform::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransform*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimHoughTransform::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransform*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimHoughTransform::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransform*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimHoughTransform::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransform*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimHoughTransformArray::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimHoughTransformArray::Class_Name()
{
   return "WCSimHoughTransformArray";
}

//______________________________________________________________________________
const char *WCSimHoughTransformArray::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransformArray*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimHoughTransformArray::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransformArray*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimHoughTransformArray::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransformArray*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimHoughTransformArray::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimHoughTransformArray*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimNtupleFactory::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimNtupleFactory::Class_Name()
{
   return "WCSimNtupleFactory";
}

//______________________________________________________________________________
const char *WCSimNtupleFactory::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleFactory*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimNtupleFactory::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleFactory*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimNtupleFactory::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleFactory*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimNtupleFactory::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleFactory*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimNtupleWriter::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimNtupleWriter::Class_Name()
{
   return "WCSimNtupleWriter";
}

//______________________________________________________________________________
const char *WCSimNtupleWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimNtupleWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimNtupleWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleWriter*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimNtupleWriter::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtupleWriter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimNtuple::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimNtuple::Class_Name()
{
   return "WCSimNtuple";
}

//______________________________________________________________________________
const char *WCSimNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtuple*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtuple*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtuple*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimNtuple::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimNtuple*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoNtuple::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoNtuple::Class_Name()
{
   return "WCSimRecoNtuple";
}

//______________________________________________________________________________
const char *WCSimRecoNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoNtuple*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoNtuple*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimRecoNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoNtuple*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimRecoNtuple::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoNtuple*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimVertexNtuple::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimVertexNtuple::Class_Name()
{
   return "WCSimVertexNtuple";
}

//______________________________________________________________________________
const char *WCSimVertexNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexNtuple*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimVertexNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexNtuple*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimVertexNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexNtuple*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimVertexNtuple::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexNtuple*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimVertexSeedNtuple::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimVertexSeedNtuple::Class_Name()
{
   return "WCSimVertexSeedNtuple";
}

//______________________________________________________________________________
const char *WCSimVertexSeedNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexSeedNtuple*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimVertexSeedNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexSeedNtuple*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimVertexSeedNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexSeedNtuple*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimVertexSeedNtuple::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimVertexSeedNtuple*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimMsg::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimMsg::Class_Name()
{
   return "WCSimMsg";
}

//______________________________________________________________________________
const char *WCSimMsg::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimMsg*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimMsg::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimMsg*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimMsg::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimMsg*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimMsg::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimMsg*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimTotalLikelihood::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimTotalLikelihood::Class_Name()
{
   return "WCSimTotalLikelihood";
}

//______________________________________________________________________________
const char *WCSimTotalLikelihood::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTotalLikelihood*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimTotalLikelihood::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTotalLikelihood*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimTotalLikelihood::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTotalLikelihood*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimTotalLikelihood::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTotalLikelihood*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimLikelihoodFitter::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimLikelihoodFitter::Class_Name()
{
   return "WCSimLikelihoodFitter";
}

//______________________________________________________________________________
const char *WCSimLikelihoodFitter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodFitter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimLikelihoodFitter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodFitter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimLikelihoodFitter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodFitter*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimLikelihoodFitter::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodFitter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimLikelihoodDigit::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimLikelihoodDigit::Class_Name()
{
   return "WCSimLikelihoodDigit";
}

//______________________________________________________________________________
const char *WCSimLikelihoodDigit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimLikelihoodDigit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimLikelihoodDigit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigit*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimLikelihoodDigit::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimLikelihoodDigitArray::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimLikelihoodDigitArray::Class_Name()
{
   return "WCSimLikelihoodDigitArray";
}

//______________________________________________________________________________
const char *WCSimLikelihoodDigitArray::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigitArray*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimLikelihoodDigitArray::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigitArray*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimLikelihoodDigitArray::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigitArray*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimLikelihoodDigitArray::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodDigitArray*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimLikelihoodTrack::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimLikelihoodTrack::Class_Name()
{
   return "WCSimLikelihoodTrack";
}

//______________________________________________________________________________
const char *WCSimLikelihoodTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimLikelihoodTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void WCSimLikelihoodTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodTrack*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *WCSimLikelihoodTrack::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void WCSimDisplay::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimDisplay.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimDisplay::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimDisplay::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimDisplay::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimDisplay::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimDisplay(void *p) {
      return  p ? new(p) ::WCSimDisplay : new ::WCSimDisplay;
   }
   static void *newArray_WCSimDisplay(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimDisplay[nElements] : new ::WCSimDisplay[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimDisplay(void *p) {
      delete ((::WCSimDisplay*)p);
   }
   static void deleteArray_WCSimDisplay(void *p) {
      delete [] ((::WCSimDisplay*)p);
   }
   static void destruct_WCSimDisplay(void *p) {
      typedef ::WCSimDisplay current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimDisplay

//______________________________________________________________________________
void WCSimDisplayFactory::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimDisplayFactory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimDisplayFactory::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimDisplayFactory::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimDisplayFactory::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimDisplayFactory::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimDisplayFactory

//______________________________________________________________________________
void WCSimDisplayViewer::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimDisplayViewer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimDisplayViewer::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimDisplayViewer::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimDisplayViewer::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimDisplayViewer::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimDisplayViewer(void *p) {
      return  p ? new(p) ::WCSimDisplayViewer : new ::WCSimDisplayViewer;
   }
   static void *newArray_WCSimDisplayViewer(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimDisplayViewer[nElements] : new ::WCSimDisplayViewer[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimDisplayViewer(void *p) {
      delete ((::WCSimDisplayViewer*)p);
   }
   static void deleteArray_WCSimDisplayViewer(void *p) {
      delete [] ((::WCSimDisplayViewer*)p);
   }
   static void destruct_WCSimDisplayViewer(void *p) {
      typedef ::WCSimDisplayViewer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimDisplayViewer

//______________________________________________________________________________
void WCSimDisplayAB::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimDisplayAB.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimDisplayAB::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimDisplayAB::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimDisplayAB::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimDisplayAB::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimDisplayAB(void *p) {
      return  p ? new(p) ::WCSimDisplayAB : new ::WCSimDisplayAB;
   }
   static void *newArray_WCSimDisplayAB(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimDisplayAB[nElements] : new ::WCSimDisplayAB[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimDisplayAB(void *p) {
      delete ((::WCSimDisplayAB*)p);
   }
   static void deleteArray_WCSimDisplayAB(void *p) {
      delete [] ((::WCSimDisplayAB*)p);
   }
   static void destruct_WCSimDisplayAB(void *p) {
      typedef ::WCSimDisplayAB current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimDisplayAB

//______________________________________________________________________________
void WCSimEveDisplay::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimEveDisplay.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimEveDisplay::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimEveDisplay::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimEveDisplay::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimEveDisplay::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimEveDisplay(void *p) {
      return  p ? new(p) ::WCSimEveDisplay : new ::WCSimEveDisplay;
   }
   static void *newArray_WCSimEveDisplay(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimEveDisplay[nElements] : new ::WCSimEveDisplay[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimEveDisplay(void *p) {
      delete ((::WCSimEveDisplay*)p);
   }
   static void deleteArray_WCSimEveDisplay(void *p) {
      delete [] ((::WCSimEveDisplay*)p);
   }
   static void destruct_WCSimEveDisplay(void *p) {
      typedef ::WCSimEveDisplay current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimEveDisplay

//______________________________________________________________________________
void WCSimVertexViewer::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimVertexViewer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimVertexViewer::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimVertexViewer::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimVertexViewer::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimVertexViewer::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimVertexViewer(void *p) {
      return  p ? new(p) ::WCSimVertexViewer : new ::WCSimVertexViewer;
   }
   static void *newArray_WCSimVertexViewer(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimVertexViewer[nElements] : new ::WCSimVertexViewer[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimVertexViewer(void *p) {
      delete ((::WCSimVertexViewer*)p);
   }
   static void deleteArray_WCSimVertexViewer(void *p) {
      delete [] ((::WCSimVertexViewer*)p);
   }
   static void destruct_WCSimVertexViewer(void *p) {
      typedef ::WCSimVertexViewer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimVertexViewer

//______________________________________________________________________________
void WCSimRingViewer::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRingViewer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRingViewer::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRingViewer::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRingViewer::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRingViewer::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRingViewer(void *p) {
      return  p ? new(p) ::WCSimRingViewer : new ::WCSimRingViewer;
   }
   static void *newArray_WCSimRingViewer(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRingViewer[nElements] : new ::WCSimRingViewer[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRingViewer(void *p) {
      delete ((::WCSimRingViewer*)p);
   }
   static void deleteArray_WCSimRingViewer(void *p) {
      delete [] ((::WCSimRingViewer*)p);
   }
   static void destruct_WCSimRingViewer(void *p) {
      typedef ::WCSimRingViewer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRingViewer

//______________________________________________________________________________
void WCSimEventWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimEventWriter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimEventWriter::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimEventWriter::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimEventWriter::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimEventWriter::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimEventWriter(void *p) {
      return  p ? new(p) ::WCSimEventWriter : new ::WCSimEventWriter;
   }
   static void *newArray_WCSimEventWriter(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimEventWriter[nElements] : new ::WCSimEventWriter[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimEventWriter(void *p) {
      delete ((::WCSimEventWriter*)p);
   }
   static void deleteArray_WCSimEventWriter(void *p) {
      delete [] ((::WCSimEventWriter*)p);
   }
   static void destruct_WCSimEventWriter(void *p) {
      typedef ::WCSimEventWriter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimEventWriter

//______________________________________________________________________________
void WCSimGeometry::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimGeometry.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimGeometry::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimGeometry::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimGeometry::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimGeometry::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimGeometry

//______________________________________________________________________________
void WCSimInterface::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimInterface.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimInterface::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimInterface::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimInterface::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimInterface::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimInterface

//______________________________________________________________________________
void WCSimParameters::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimParameters.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimParameters::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimParameters::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimParameters::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimParameters::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimParameters

//______________________________________________________________________________
void WCSimTrueEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimTrueEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimTrueEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimTrueEvent::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimTrueEvent::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimTrueEvent::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimTrueEvent(void *p) {
      return  p ? new(p) ::WCSimTrueEvent : new ::WCSimTrueEvent;
   }
   static void *newArray_WCSimTrueEvent(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimTrueEvent[nElements] : new ::WCSimTrueEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimTrueEvent(void *p) {
      delete ((::WCSimTrueEvent*)p);
   }
   static void deleteArray_WCSimTrueEvent(void *p) {
      delete [] ((::WCSimTrueEvent*)p);
   }
   static void destruct_WCSimTrueEvent(void *p) {
      typedef ::WCSimTrueEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTrueEvent

//______________________________________________________________________________
void WCSimTrueTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimTrueTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimTrueTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimTrueTrack::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimTrueTrack::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimTrueTrack::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimTrueTrack(void *p) {
      delete ((::WCSimTrueTrack*)p);
   }
   static void deleteArray_WCSimTrueTrack(void *p) {
      delete [] ((::WCSimTrueTrack*)p);
   }
   static void destruct_WCSimTrueTrack(void *p) {
      typedef ::WCSimTrueTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTrueTrack

//______________________________________________________________________________
void WCSimRecoObjectTable::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoObjectTable.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoObjectTable::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoObjectTable::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoObjectTable::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoObjectTable::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimRecoObjectTable

//______________________________________________________________________________
void WCSimRecoFactory::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoFactory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoFactory::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoFactory::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoFactory::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoFactory::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimRecoFactory

//______________________________________________________________________________
void WCSimReco::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimReco.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimReco::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimReco::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimReco::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimReco::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimReco(void *p) {
      delete ((::WCSimReco*)p);
   }
   static void deleteArray_WCSimReco(void *p) {
      delete [] ((::WCSimReco*)p);
   }
   static void destruct_WCSimReco(void *p) {
      typedef ::WCSimReco current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimReco

//______________________________________________________________________________
void WCSimRecoAB::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoAB.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoAB::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoAB::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoAB::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoAB::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRecoAB(void *p) {
      return  p ? new(p) ::WCSimRecoAB : new ::WCSimRecoAB;
   }
   static void *newArray_WCSimRecoAB(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRecoAB[nElements] : new ::WCSimRecoAB[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRecoAB(void *p) {
      delete ((::WCSimRecoAB*)p);
   }
   static void deleteArray_WCSimRecoAB(void *p) {
      delete [] ((::WCSimRecoAB*)p);
   }
   static void destruct_WCSimRecoAB(void *p) {
      typedef ::WCSimRecoAB current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoAB

//______________________________________________________________________________
void WCSimRecoEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoEvent::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoEvent::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoEvent::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRecoEvent(void *p) {
      return  p ? new(p) ::WCSimRecoEvent : new ::WCSimRecoEvent;
   }
   static void *newArray_WCSimRecoEvent(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRecoEvent[nElements] : new ::WCSimRecoEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRecoEvent(void *p) {
      delete ((::WCSimRecoEvent*)p);
   }
   static void deleteArray_WCSimRecoEvent(void *p) {
      delete [] ((::WCSimRecoEvent*)p);
   }
   static void destruct_WCSimRecoEvent(void *p) {
      typedef ::WCSimRecoEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoEvent

//______________________________________________________________________________
void WCSimRecoDigit::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoDigit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoDigit::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoDigit::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoDigit::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoDigit::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimRecoDigit(void *p) {
      delete ((::WCSimRecoDigit*)p);
   }
   static void deleteArray_WCSimRecoDigit(void *p) {
      delete [] ((::WCSimRecoDigit*)p);
   }
   static void destruct_WCSimRecoDigit(void *p) {
      typedef ::WCSimRecoDigit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoDigit

//______________________________________________________________________________
void WCSimRecoCluster::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoCluster.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoCluster::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoCluster::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoCluster::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoCluster::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRecoCluster(void *p) {
      return  p ? new(p) ::WCSimRecoCluster : new ::WCSimRecoCluster;
   }
   static void *newArray_WCSimRecoCluster(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRecoCluster[nElements] : new ::WCSimRecoCluster[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRecoCluster(void *p) {
      delete ((::WCSimRecoCluster*)p);
   }
   static void deleteArray_WCSimRecoCluster(void *p) {
      delete [] ((::WCSimRecoCluster*)p);
   }
   static void destruct_WCSimRecoCluster(void *p) {
      typedef ::WCSimRecoCluster current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoCluster

//______________________________________________________________________________
void WCSimRecoClusterDigit::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoClusterDigit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoClusterDigit::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoClusterDigit::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoClusterDigit::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoClusterDigit::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimRecoClusterDigit(void *p) {
      delete ((::WCSimRecoClusterDigit*)p);
   }
   static void deleteArray_WCSimRecoClusterDigit(void *p) {
      delete [] ((::WCSimRecoClusterDigit*)p);
   }
   static void destruct_WCSimRecoClusterDigit(void *p) {
      typedef ::WCSimRecoClusterDigit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoClusterDigit

//______________________________________________________________________________
void WCSimRecoRing::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoRing.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoRing::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoRing::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoRing::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoRing::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimRecoRing(void *p) {
      delete ((::WCSimRecoRing*)p);
   }
   static void deleteArray_WCSimRecoRing(void *p) {
      delete [] ((::WCSimRecoRing*)p);
   }
   static void destruct_WCSimRecoRing(void *p) {
      typedef ::WCSimRecoRing current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoRing

//______________________________________________________________________________
void WCSimRecoVertex::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoVertex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoVertex::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoVertex::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoVertex::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoVertex::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRecoVertex(void *p) {
      return  p ? new(p) ::WCSimRecoVertex : new ::WCSimRecoVertex;
   }
   static void *newArray_WCSimRecoVertex(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRecoVertex[nElements] : new ::WCSimRecoVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRecoVertex(void *p) {
      delete ((::WCSimRecoVertex*)p);
   }
   static void deleteArray_WCSimRecoVertex(void *p) {
      delete [] ((::WCSimRecoVertex*)p);
   }
   static void destruct_WCSimRecoVertex(void *p) {
      typedef ::WCSimRecoVertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoVertex

//______________________________________________________________________________
void WCSimVertexFinder::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimVertexFinder.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimVertexFinder::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimVertexFinder::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimVertexFinder::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimVertexFinder::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimVertexFinder

//______________________________________________________________________________
void WCSimRingFinder::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRingFinder.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRingFinder::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRingFinder::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRingFinder::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRingFinder::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimRingFinder

//______________________________________________________________________________
void WCSimDataCleaner::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimDataCleaner.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimDataCleaner::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimDataCleaner::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimDataCleaner::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimDataCleaner::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimDataCleaner

//______________________________________________________________________________
void WCSimVertexGeometry::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimVertexGeometry.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimVertexGeometry::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimVertexGeometry::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimVertexGeometry::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimVertexGeometry::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimVertexGeometry

//______________________________________________________________________________
void WCSimHoughTransform::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimHoughTransform.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimHoughTransform::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimHoughTransform::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimHoughTransform::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimHoughTransform::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimHoughTransform(void *p) {
      delete ((::WCSimHoughTransform*)p);
   }
   static void deleteArray_WCSimHoughTransform(void *p) {
      delete [] ((::WCSimHoughTransform*)p);
   }
   static void destruct_WCSimHoughTransform(void *p) {
      typedef ::WCSimHoughTransform current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimHoughTransform

//______________________________________________________________________________
void WCSimHoughTransformArray::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimHoughTransformArray.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimHoughTransformArray::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimHoughTransformArray::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimHoughTransformArray::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimHoughTransformArray::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimHoughTransformArray(void *p) {
      delete ((::WCSimHoughTransformArray*)p);
   }
   static void deleteArray_WCSimHoughTransformArray(void *p) {
      delete [] ((::WCSimHoughTransformArray*)p);
   }
   static void destruct_WCSimHoughTransformArray(void *p) {
      typedef ::WCSimHoughTransformArray current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimHoughTransformArray

//______________________________________________________________________________
void WCSimNtupleFactory::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimNtupleFactory.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimNtupleFactory::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimNtupleFactory::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimNtupleFactory::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimNtupleFactory::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimNtupleFactory

//______________________________________________________________________________
void WCSimNtupleWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimNtupleWriter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimNtupleWriter::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimNtupleWriter::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimNtupleWriter::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimNtupleWriter::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimNtupleWriter(void *p) {
      return  p ? new(p) ::WCSimNtupleWriter : new ::WCSimNtupleWriter;
   }
   static void *newArray_WCSimNtupleWriter(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimNtupleWriter[nElements] : new ::WCSimNtupleWriter[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimNtupleWriter(void *p) {
      delete ((::WCSimNtupleWriter*)p);
   }
   static void deleteArray_WCSimNtupleWriter(void *p) {
      delete [] ((::WCSimNtupleWriter*)p);
   }
   static void destruct_WCSimNtupleWriter(void *p) {
      typedef ::WCSimNtupleWriter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimNtupleWriter

//______________________________________________________________________________
void WCSimNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimNtuple.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimNtuple::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimNtuple::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimNtuple::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimNtuple::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimNtuple(void *p) {
      return  p ? new(p) ::WCSimNtuple : new ::WCSimNtuple;
   }
   static void *newArray_WCSimNtuple(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimNtuple[nElements] : new ::WCSimNtuple[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimNtuple(void *p) {
      delete ((::WCSimNtuple*)p);
   }
   static void deleteArray_WCSimNtuple(void *p) {
      delete [] ((::WCSimNtuple*)p);
   }
   static void destruct_WCSimNtuple(void *p) {
      typedef ::WCSimNtuple current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimNtuple

//______________________________________________________________________________
void WCSimRecoNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoNtuple.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoNtuple::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoNtuple::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimRecoNtuple::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimRecoNtuple::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRecoNtuple(void *p) {
      return  p ? new(p) ::WCSimRecoNtuple : new ::WCSimRecoNtuple;
   }
   static void *newArray_WCSimRecoNtuple(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRecoNtuple[nElements] : new ::WCSimRecoNtuple[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRecoNtuple(void *p) {
      delete ((::WCSimRecoNtuple*)p);
   }
   static void deleteArray_WCSimRecoNtuple(void *p) {
      delete [] ((::WCSimRecoNtuple*)p);
   }
   static void destruct_WCSimRecoNtuple(void *p) {
      typedef ::WCSimRecoNtuple current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoNtuple

//______________________________________________________________________________
void WCSimVertexNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimVertexNtuple.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimVertexNtuple::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimVertexNtuple::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimVertexNtuple::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimVertexNtuple::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimVertexNtuple(void *p) {
      return  p ? new(p) ::WCSimVertexNtuple : new ::WCSimVertexNtuple;
   }
   static void *newArray_WCSimVertexNtuple(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimVertexNtuple[nElements] : new ::WCSimVertexNtuple[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimVertexNtuple(void *p) {
      delete ((::WCSimVertexNtuple*)p);
   }
   static void deleteArray_WCSimVertexNtuple(void *p) {
      delete [] ((::WCSimVertexNtuple*)p);
   }
   static void destruct_WCSimVertexNtuple(void *p) {
      typedef ::WCSimVertexNtuple current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimVertexNtuple

//______________________________________________________________________________
void WCSimVertexSeedNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimVertexSeedNtuple.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimVertexSeedNtuple::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimVertexSeedNtuple::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimVertexSeedNtuple::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimVertexSeedNtuple::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimVertexSeedNtuple(void *p) {
      return  p ? new(p) ::WCSimVertexSeedNtuple : new ::WCSimVertexSeedNtuple;
   }
   static void *newArray_WCSimVertexSeedNtuple(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimVertexSeedNtuple[nElements] : new ::WCSimVertexSeedNtuple[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimVertexSeedNtuple(void *p) {
      delete ((::WCSimVertexSeedNtuple*)p);
   }
   static void deleteArray_WCSimVertexSeedNtuple(void *p) {
      delete [] ((::WCSimVertexSeedNtuple*)p);
   }
   static void destruct_WCSimVertexSeedNtuple(void *p) {
      typedef ::WCSimVertexSeedNtuple current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimVertexSeedNtuple

//______________________________________________________________________________
void WCSimMsg::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimMsg.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimMsg::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimMsg::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimMsg::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimMsg::IsA());
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimMsg

//______________________________________________________________________________
void WCSimTotalLikelihood::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimTotalLikelihood.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimTotalLikelihood::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimTotalLikelihood::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimTotalLikelihood::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimTotalLikelihood::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimTotalLikelihood(void *p) {
      delete ((::WCSimTotalLikelihood*)p);
   }
   static void deleteArray_WCSimTotalLikelihood(void *p) {
      delete [] ((::WCSimTotalLikelihood*)p);
   }
   static void destruct_WCSimTotalLikelihood(void *p) {
      typedef ::WCSimTotalLikelihood current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTotalLikelihood

//______________________________________________________________________________
void WCSimLikelihoodFitter::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimLikelihoodFitter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimLikelihoodFitter::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimLikelihoodFitter::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimLikelihoodFitter::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimLikelihoodFitter::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimLikelihoodFitter(void *p) {
      delete ((::WCSimLikelihoodFitter*)p);
   }
   static void deleteArray_WCSimLikelihoodFitter(void *p) {
      delete [] ((::WCSimLikelihoodFitter*)p);
   }
   static void destruct_WCSimLikelihoodFitter(void *p) {
      typedef ::WCSimLikelihoodFitter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimLikelihoodFitter

//______________________________________________________________________________
namespace ROOT {
   void WCSimLikelihoodTuner_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
   gInterpreter->InspectMembers(R__insp, obj, ::ROOT::GenerateInitInstanceLocal((const ::WCSimLikelihoodTuner*)0x0)->GetClass());
   }

}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimLikelihoodTuner(void *p) {
      delete ((::WCSimLikelihoodTuner*)p);
   }
   static void deleteArray_WCSimLikelihoodTuner(void *p) {
      delete [] ((::WCSimLikelihoodTuner*)p);
   }
   static void destruct_WCSimLikelihoodTuner(void *p) {
      typedef ::WCSimLikelihoodTuner current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimLikelihoodTuner

//______________________________________________________________________________
namespace ROOT {
   void WCSimChargeLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
   gInterpreter->InspectMembers(R__insp, obj, ::ROOT::GenerateInitInstanceLocal((const ::WCSimChargeLikelihood*)0x0)->GetClass());
   }

}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimChargeLikelihood(void *p) {
      delete ((::WCSimChargeLikelihood*)p);
   }
   static void deleteArray_WCSimChargeLikelihood(void *p) {
      delete [] ((::WCSimChargeLikelihood*)p);
   }
   static void destruct_WCSimChargeLikelihood(void *p) {
      typedef ::WCSimChargeLikelihood current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimChargeLikelihood

//______________________________________________________________________________
void WCSimLikelihoodDigit::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimLikelihoodDigit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimLikelihoodDigit::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimLikelihoodDigit::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimLikelihoodDigit::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimLikelihoodDigit::IsA());
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimLikelihoodDigit(void *p) {
      delete ((::WCSimLikelihoodDigit*)p);
   }
   static void deleteArray_WCSimLikelihoodDigit(void *p) {
      delete [] ((::WCSimLikelihoodDigit*)p);
   }
   static void destruct_WCSimLikelihoodDigit(void *p) {
      typedef ::WCSimLikelihoodDigit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimLikelihoodDigit

//______________________________________________________________________________
void WCSimLikelihoodDigitArray::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimLikelihoodDigitArray.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimLikelihoodDigitArray::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimLikelihoodDigitArray::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimLikelihoodDigitArray::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimLikelihoodDigitArray::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimLikelihoodDigitArray(void *p) {
      return  p ? new(p) ::WCSimLikelihoodDigitArray : new ::WCSimLikelihoodDigitArray;
   }
   static void *newArray_WCSimLikelihoodDigitArray(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimLikelihoodDigitArray[nElements] : new ::WCSimLikelihoodDigitArray[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimLikelihoodDigitArray(void *p) {
      delete ((::WCSimLikelihoodDigitArray*)p);
   }
   static void deleteArray_WCSimLikelihoodDigitArray(void *p) {
      delete [] ((::WCSimLikelihoodDigitArray*)p);
   }
   static void destruct_WCSimLikelihoodDigitArray(void *p) {
      typedef ::WCSimLikelihoodDigitArray current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimLikelihoodDigitArray

//______________________________________________________________________________
void WCSimLikelihoodTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimLikelihoodTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimLikelihoodTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimLikelihoodTrack::Class(),this);
   }
}

//______________________________________________________________________________
void WCSimLikelihoodTrack::ShowMembers(TMemberInspector &R__insp)
{
   gInterpreter->InspectMembers(R__insp, this, ::WCSimLikelihoodTrack::IsA());
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimLikelihoodTrack(void *p) {
      return  p ? new(p) ::WCSimLikelihoodTrack : new ::WCSimLikelihoodTrack;
   }
   static void *newArray_WCSimLikelihoodTrack(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimLikelihoodTrack[nElements] : new ::WCSimLikelihoodTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimLikelihoodTrack(void *p) {
      delete ((::WCSimLikelihoodTrack*)p);
   }
   static void deleteArray_WCSimLikelihoodTrack(void *p) {
      delete [] ((::WCSimLikelihoodTrack*)p);
   }
   static void destruct_WCSimLikelihoodTrack(void *p) {
      typedef ::WCSimLikelihoodTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimLikelihoodTrack

//______________________________________________________________________________
namespace ROOT {
   void WCSimTimeLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
   gInterpreter->InspectMembers(R__insp, obj, ::ROOT::GenerateInitInstanceLocal((const ::WCSimTimeLikelihood*)0x0)->GetClass());
   }

}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimTimeLikelihood(void *p) {
      delete ((::WCSimTimeLikelihood*)p);
   }
   static void deleteArray_WCSimTimeLikelihood(void *p) {
      delete [] ((::WCSimTimeLikelihood*)p);
   }
   static void destruct_WCSimTimeLikelihood(void *p) {
      typedef ::WCSimTimeLikelihood current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTimeLikelihood

//______________________________________________________________________________
namespace ROOT {
   void WCSimDigitizerLikelihood_ShowMembers(void *obj, TMemberInspector &R__insp)
   {
   gInterpreter->InspectMembers(R__insp, obj, ::ROOT::GenerateInitInstanceLocal((const ::WCSimDigitizerLikelihood*)0x0)->GetClass());
   }

}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimDigitizerLikelihood(void *p) {
      return  p ? new(p) ::WCSimDigitizerLikelihood : new ::WCSimDigitizerLikelihood;
   }
   static void *newArray_WCSimDigitizerLikelihood(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimDigitizerLikelihood[nElements] : new ::WCSimDigitizerLikelihood[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimDigitizerLikelihood(void *p) {
      delete ((::WCSimDigitizerLikelihood*)p);
   }
   static void deleteArray_WCSimDigitizerLikelihood(void *p) {
      delete [] ((::WCSimDigitizerLikelihood*)p);
   }
   static void destruct_WCSimDigitizerLikelihood(void *p) {
      typedef ::WCSimDigitizerLikelihood current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimDigitizerLikelihood

namespace ROOT {
   void vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_Dictionary();
   static void vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p);
   static void destruct_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<WCSimLikelihoodTrack> >*)
   {
      vector<vector<WCSimLikelihoodTrack> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<WCSimLikelihoodTrack> >),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<WCSimLikelihoodTrack> >", -2, "vector", 208,
                  typeid(vector<vector<WCSimLikelihoodTrack> >), DefineBehavior(ptr, ptr),
                  0, &vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<WCSimLikelihoodTrack> >) );
      instance.SetNew(&new_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<WCSimLikelihoodTrack> > >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<vector<WCSimLikelihoodTrack> >*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<WCSimLikelihoodTrack> >*)0x0)->GetClass();
      vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_TClassManip(theClass);
   }

   static void vectorlEvectorlEWCSimLikelihoodTrackgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<vector<WCSimLikelihoodTrack> > : new vector<vector<WCSimLikelihoodTrack> >;
   }
   static void *newArray_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<vector<WCSimLikelihoodTrack> >[nElements] : new vector<vector<WCSimLikelihoodTrack> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p) {
      delete ((vector<vector<WCSimLikelihoodTrack> >*)p);
   }
   static void deleteArray_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p) {
      delete [] ((vector<vector<WCSimLikelihoodTrack> >*)p);
   }
   static void destruct_vectorlEvectorlEWCSimLikelihoodTrackgRsPgR(void *p) {
      typedef vector<vector<WCSimLikelihoodTrack> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<WCSimLikelihoodTrack> >

namespace ROOT {
   void vectorlEdoublegR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 208,
                  typeid(vector<double>), DefineBehavior(ptr, ptr),
                  0, &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   void vectorlEWCSimRecoDigitmUgR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEWCSimRecoDigitmUgR_Dictionary();
   static void vectorlEWCSimRecoDigitmUgR_TClassManip(TClass*);
   static void *new_vectorlEWCSimRecoDigitmUgR(void *p = 0);
   static void *newArray_vectorlEWCSimRecoDigitmUgR(Long_t size, void *p);
   static void delete_vectorlEWCSimRecoDigitmUgR(void *p);
   static void deleteArray_vectorlEWCSimRecoDigitmUgR(void *p);
   static void destruct_vectorlEWCSimRecoDigitmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<WCSimRecoDigit*>*)
   {
      vector<WCSimRecoDigit*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<WCSimRecoDigit*>),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<WCSimRecoDigit*>", -2, "vector", 208,
                  typeid(vector<WCSimRecoDigit*>), DefineBehavior(ptr, ptr),
                  0, &vectorlEWCSimRecoDigitmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<WCSimRecoDigit*>) );
      instance.SetNew(&new_vectorlEWCSimRecoDigitmUgR);
      instance.SetNewArray(&newArray_vectorlEWCSimRecoDigitmUgR);
      instance.SetDelete(&delete_vectorlEWCSimRecoDigitmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEWCSimRecoDigitmUgR);
      instance.SetDestructor(&destruct_vectorlEWCSimRecoDigitmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<WCSimRecoDigit*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<WCSimRecoDigit*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEWCSimRecoDigitmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<WCSimRecoDigit*>*)0x0)->GetClass();
      vectorlEWCSimRecoDigitmUgR_TClassManip(theClass);
   }

   static void vectorlEWCSimRecoDigitmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEWCSimRecoDigitmUgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<WCSimRecoDigit*> : new vector<WCSimRecoDigit*>;
   }
   static void *newArray_vectorlEWCSimRecoDigitmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<WCSimRecoDigit*>[nElements] : new vector<WCSimRecoDigit*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEWCSimRecoDigitmUgR(void *p) {
      delete ((vector<WCSimRecoDigit*>*)p);
   }
   static void deleteArray_vectorlEWCSimRecoDigitmUgR(void *p) {
      delete [] ((vector<WCSimRecoDigit*>*)p);
   }
   static void destruct_vectorlEWCSimRecoDigitmUgR(void *p) {
      typedef vector<WCSimRecoDigit*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<WCSimRecoDigit*>

namespace ROOT {
   void vectorlEWCSimLikelihoodTrackgR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEWCSimLikelihoodTrackgR_Dictionary();
   static void vectorlEWCSimLikelihoodTrackgR_TClassManip(TClass*);
   static void *new_vectorlEWCSimLikelihoodTrackgR(void *p = 0);
   static void *newArray_vectorlEWCSimLikelihoodTrackgR(Long_t size, void *p);
   static void delete_vectorlEWCSimLikelihoodTrackgR(void *p);
   static void deleteArray_vectorlEWCSimLikelihoodTrackgR(void *p);
   static void destruct_vectorlEWCSimLikelihoodTrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<WCSimLikelihoodTrack>*)
   {
      vector<WCSimLikelihoodTrack> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<WCSimLikelihoodTrack>),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<WCSimLikelihoodTrack>", -2, "vector", 208,
                  typeid(vector<WCSimLikelihoodTrack>), DefineBehavior(ptr, ptr),
                  0, &vectorlEWCSimLikelihoodTrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<WCSimLikelihoodTrack>) );
      instance.SetNew(&new_vectorlEWCSimLikelihoodTrackgR);
      instance.SetNewArray(&newArray_vectorlEWCSimLikelihoodTrackgR);
      instance.SetDelete(&delete_vectorlEWCSimLikelihoodTrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEWCSimLikelihoodTrackgR);
      instance.SetDestructor(&destruct_vectorlEWCSimLikelihoodTrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<WCSimLikelihoodTrack> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<WCSimLikelihoodTrack>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEWCSimLikelihoodTrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<WCSimLikelihoodTrack>*)0x0)->GetClass();
      vectorlEWCSimLikelihoodTrackgR_TClassManip(theClass);
   }

   static void vectorlEWCSimLikelihoodTrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEWCSimLikelihoodTrackgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<WCSimLikelihoodTrack> : new vector<WCSimLikelihoodTrack>;
   }
   static void *newArray_vectorlEWCSimLikelihoodTrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<WCSimLikelihoodTrack>[nElements] : new vector<WCSimLikelihoodTrack>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEWCSimLikelihoodTrackgR(void *p) {
      delete ((vector<WCSimLikelihoodTrack>*)p);
   }
   static void deleteArray_vectorlEWCSimLikelihoodTrackgR(void *p) {
      delete [] ((vector<WCSimLikelihoodTrack>*)p);
   }
   static void destruct_vectorlEWCSimLikelihoodTrackgR(void *p) {
      typedef vector<WCSimLikelihoodTrack> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<WCSimLikelihoodTrack>

namespace ROOT {
   void vectorlEWCSimLikelihoodTrackmUgR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void vectorlEWCSimLikelihoodTrackmUgR_Dictionary();
   static void vectorlEWCSimLikelihoodTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlEWCSimLikelihoodTrackmUgR(void *p = 0);
   static void *newArray_vectorlEWCSimLikelihoodTrackmUgR(Long_t size, void *p);
   static void delete_vectorlEWCSimLikelihoodTrackmUgR(void *p);
   static void deleteArray_vectorlEWCSimLikelihoodTrackmUgR(void *p);
   static void destruct_vectorlEWCSimLikelihoodTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<WCSimLikelihoodTrack*>*)
   {
      vector<WCSimLikelihoodTrack*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<WCSimLikelihoodTrack*>),0);
      static ::ROOT::TGenericClassInfo 
         instance("vector<WCSimLikelihoodTrack*>", -2, "vector", 208,
                  typeid(vector<WCSimLikelihoodTrack*>), DefineBehavior(ptr, ptr),
                  0, &vectorlEWCSimLikelihoodTrackmUgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<WCSimLikelihoodTrack*>) );
      instance.SetNew(&new_vectorlEWCSimLikelihoodTrackmUgR);
      instance.SetNewArray(&newArray_vectorlEWCSimLikelihoodTrackmUgR);
      instance.SetDelete(&delete_vectorlEWCSimLikelihoodTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEWCSimLikelihoodTrackmUgR);
      instance.SetDestructor(&destruct_vectorlEWCSimLikelihoodTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<WCSimLikelihoodTrack*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<WCSimLikelihoodTrack*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void vectorlEWCSimLikelihoodTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<WCSimLikelihoodTrack*>*)0x0)->GetClass();
      vectorlEWCSimLikelihoodTrackmUgR_TClassManip(theClass);
   }

   static void vectorlEWCSimLikelihoodTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEWCSimLikelihoodTrackmUgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<WCSimLikelihoodTrack*> : new vector<WCSimLikelihoodTrack*>;
   }
   static void *newArray_vectorlEWCSimLikelihoodTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<WCSimLikelihoodTrack*>[nElements] : new vector<WCSimLikelihoodTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEWCSimLikelihoodTrackmUgR(void *p) {
      delete ((vector<WCSimLikelihoodTrack*>*)p);
   }
   static void deleteArray_vectorlEWCSimLikelihoodTrackmUgR(void *p) {
      delete [] ((vector<WCSimLikelihoodTrack*>*)p);
   }
   static void destruct_vectorlEWCSimLikelihoodTrackmUgR(void *p) {
      typedef vector<WCSimLikelihoodTrack*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<WCSimLikelihoodTrack*>

namespace ROOT {
   void maplEintcOintgR_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void maplEintcOintgR_Dictionary();
   static void maplEintcOintgR_TClassManip(TClass*);
   static void *new_maplEintcOintgR(void *p = 0);
   static void *newArray_maplEintcOintgR(Long_t size, void *p);
   static void delete_maplEintcOintgR(void *p);
   static void deleteArray_maplEintcOintgR(void *p);
   static void destruct_maplEintcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,int>*)
   {
      map<int,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,int>),0);
      static ::ROOT::TGenericClassInfo 
         instance("map<int,int>", -2, "map", 90,
                  typeid(map<int,int>), DefineBehavior(ptr, ptr),
                  0, &maplEintcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,int>) );
      instance.SetNew(&new_maplEintcOintgR);
      instance.SetNewArray(&newArray_maplEintcOintgR);
      instance.SetDelete(&delete_maplEintcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOintgR);
      instance.SetDestructor(&destruct_maplEintcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<int,int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static void maplEintcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,int>*)0x0)->GetClass();
      maplEintcOintgR_TClassManip(theClass);
   }

   static void maplEintcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOintgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) map<int,int> : new map<int,int>;
   }
   static void *newArray_maplEintcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) map<int,int>[nElements] : new map<int,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOintgR(void *p) {
      delete ((map<int,int>*)p);
   }
   static void deleteArray_maplEintcOintgR(void *p) {
      delete [] ((map<int,int>*)p);
   }
   static void destruct_maplEintcOintgR(void *p) {
      typedef map<int,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,int>

namespace {
  void TriggerDictionaryInitialization_WCSimAnalysisRootDict_Impl() {
    static const char* headers[] = {
"WCSimDigitizerLikelihood.hh",
"WCSimTotalLikelihood.hh",
"WCSimLikelihoodTrack.hh",
"WCSimLikelihoodDigit.hh",
"WCSimLikelihoodDigitArray.hh",
"WCSimChargeLikelihood.hh",
"WCSimLikelihoodTuner.hh",
"WCSimLikelihoodFitter.hh",
"WCSimDisplayViewer.hh",
"WCSimDisplayFactory.hh",
"WCSimDisplay.hh",
"WCSimDisplayAB.hh",
"WCSimEveDisplay.hh",
"WCSimEventWriter.hh",
"WCSimGeometry.hh",
"WCSimInterface.hh",
"WCSimParameters.hh",
"WCSimRecoObjectTable.hh",
"WCSimRecoFactory.hh",
"WCSimReco.hh",
"WCSimRecoAB.hh",
"WCSimRecoDigit.hh",
"WCSimRecoCluster.hh",
"WCSimRecoClusterDigit.hh",
"WCSimRecoRing.hh",
"WCSimRecoVertex.hh",
"WCSimRecoEvent.hh",
"WCSimTrueEvent.hh",
"WCSimTrueTrack.hh",
"WCSimHoughTransform.hh",
"WCSimHoughTransformArray.hh",
"WCSimDataCleaner.hh",
"WCSimVertexFinder.hh",
"WCSimVertexGeometry.hh",
"WCSimVertexViewer.hh",
"WCSimRingFinder.hh",
"WCSimRingViewer.hh",
"WCSimNtupleFactory.hh",
"WCSimNtuple.hh",
"WCSimRecoNtuple.hh",
"WCSimVertexNtuple.hh",
"WCSimVertexSeedNtuple.hh",
"WCSimNtupleWriter.hh",
"WCSimMsg.hh",
0
    };
    static const char* allHeaders[] = {
"././src/WCSimAnalysisRootDict7932b7084f_dictContent.h",
"././src/WCSimAnalysisRootDict7932b7084f_dictUmbrella.h",
"./include/WCSimAnalysisRootLinkDef.hh",
"./include/WCSimChargeLikelihood.hh",
"./include/WCSimDataCleaner.hh",
"./include/WCSimDigitizerLikelihood.hh",
"./include/WCSimDisplay.hh",
"./include/WCSimDisplayAB.hh",
"./include/WCSimDisplayFactory.hh",
"./include/WCSimDisplayViewer.hh",
"./include/WCSimEveDisplay.hh",
"./include/WCSimEventWriter.hh",
"./include/WCSimGeometry.hh",
"./include/WCSimHoughTransform.hh",
"./include/WCSimHoughTransformArray.hh",
"./include/WCSimInterface.hh",
"./include/WCSimLikelihoodDigit.hh",
"./include/WCSimLikelihoodDigitArray.hh",
"./include/WCSimLikelihoodFitter.hh",
"./include/WCSimLikelihoodTrack.hh",
"./include/WCSimLikelihoodTuner.hh",
"./include/WCSimMsg.hh",
"./include/WCSimNtuple.hh",
"./include/WCSimNtupleFactory.hh",
"./include/WCSimNtupleWriter.hh",
"./include/WCSimParameters.hh",
"./include/WCSimReco.hh",
"./include/WCSimRecoAB.hh",
"./include/WCSimRecoCluster.hh",
"./include/WCSimRecoClusterDigit.hh",
"./include/WCSimRecoDigit.hh",
"./include/WCSimRecoEvent.hh",
"./include/WCSimRecoFactory.hh",
"./include/WCSimRecoNtuple.hh",
"./include/WCSimRecoObjectTable.hh",
"./include/WCSimRecoRing.hh",
"./include/WCSimRecoVertex.hh",
"./include/WCSimRingFinder.hh",
"./include/WCSimRingViewer.hh",
"./include/WCSimTimeLikelihood.hh",
"./include/WCSimTotalLikelihood.hh",
"./include/WCSimTrueEvent.hh",
"./include/WCSimTrueTrack.hh",
"./include/WCSimVertexFinder.hh",
"./include/WCSimVertexGeometry.hh",
"./include/WCSimVertexNtuple.hh",
"./include/WCSimVertexSeedNtuple.hh",
"./include/WCSimVertexViewer.hh",
"/etc/root/cling/Interpreter/RuntimeException.h",
"/etc/root/cling/Interpreter/RuntimeUniverse.h",
"/etc/root/cling/Interpreter/ValuePrinter.h",
"/etc/root/cling/Interpreter/ValuePrinterInfo.h",
"/etc/root/cling/lib/clang/3.4/include/float.h",
"/etc/root/cling/lib/clang/3.4/include/stdarg.h",
"/etc/root/cling/lib/clang/3.4/include/stddef.h",
"/home/andy/work/CHIPS/WCSim/include/WCSimRootEvent.hh",
"/home/andy/work/CHIPS/WCSim/include/WCSimRootGeom.hh",
"/usr/include/_G_config.h",
"/usr/include/alloca.h",
"/usr/include/assert.h",
"/usr/include/c++/4.7/algorithm",
"/usr/include/c++/4.7/backward/binders.h",
"/usr/include/c++/4.7/bits/algorithmfwd.h",
"/usr/include/c++/4.7/bits/allocator.h",
"/usr/include/c++/4.7/bits/atomic_lockfree_defines.h",
"/usr/include/c++/4.7/bits/basic_ios.h",
"/usr/include/c++/4.7/bits/basic_ios.tcc",
"/usr/include/c++/4.7/bits/basic_string.h",
"/usr/include/c++/4.7/bits/basic_string.tcc",
"/usr/include/c++/4.7/bits/char_traits.h",
"/usr/include/c++/4.7/bits/concept_check.h",
"/usr/include/c++/4.7/bits/cpp_type_traits.h",
"/usr/include/c++/4.7/bits/cxxabi_forced.h",
"/usr/include/c++/4.7/bits/exception_defines.h",
"/usr/include/c++/4.7/bits/functexcept.h",
"/usr/include/c++/4.7/bits/ios_base.h",
"/usr/include/c++/4.7/bits/istream.tcc",
"/usr/include/c++/4.7/bits/locale_classes.h",
"/usr/include/c++/4.7/bits/locale_classes.tcc",
"/usr/include/c++/4.7/bits/locale_facets.h",
"/usr/include/c++/4.7/bits/locale_facets.tcc",
"/usr/include/c++/4.7/bits/localefwd.h",
"/usr/include/c++/4.7/bits/move.h",
"/usr/include/c++/4.7/bits/ostream.tcc",
"/usr/include/c++/4.7/bits/ostream_insert.h",
"/usr/include/c++/4.7/bits/postypes.h",
"/usr/include/c++/4.7/bits/range_access.h",
"/usr/include/c++/4.7/bits/stl_algo.h",
"/usr/include/c++/4.7/bits/stl_algobase.h",
"/usr/include/c++/4.7/bits/stl_bvector.h",
"/usr/include/c++/4.7/bits/stl_construct.h",
"/usr/include/c++/4.7/bits/stl_function.h",
"/usr/include/c++/4.7/bits/stl_heap.h",
"/usr/include/c++/4.7/bits/stl_iterator.h",
"/usr/include/c++/4.7/bits/stl_iterator_base_funcs.h",
"/usr/include/c++/4.7/bits/stl_iterator_base_types.h",
"/usr/include/c++/4.7/bits/stl_map.h",
"/usr/include/c++/4.7/bits/stl_multimap.h",
"/usr/include/c++/4.7/bits/stl_multiset.h",
"/usr/include/c++/4.7/bits/stl_pair.h",
"/usr/include/c++/4.7/bits/stl_relops.h",
"/usr/include/c++/4.7/bits/stl_set.h",
"/usr/include/c++/4.7/bits/stl_tempbuf.h",
"/usr/include/c++/4.7/bits/stl_tree.h",
"/usr/include/c++/4.7/bits/stl_uninitialized.h",
"/usr/include/c++/4.7/bits/stl_vector.h",
"/usr/include/c++/4.7/bits/stream_iterator.h",
"/usr/include/c++/4.7/bits/streambuf.tcc",
"/usr/include/c++/4.7/bits/streambuf_iterator.h",
"/usr/include/c++/4.7/bits/stringfwd.h",
"/usr/include/c++/4.7/bits/vector.tcc",
"/usr/include/c++/4.7/cctype",
"/usr/include/c++/4.7/clocale",
"/usr/include/c++/4.7/cmath",
"/usr/include/c++/4.7/cstdlib",
"/usr/include/c++/4.7/cwchar",
"/usr/include/c++/4.7/cwctype",
"/usr/include/c++/4.7/debug/debug.h",
"/usr/include/c++/4.7/exception",
"/usr/include/c++/4.7/ext/alloc_traits.h",
"/usr/include/c++/4.7/ext/atomicity.h",
"/usr/include/c++/4.7/ext/new_allocator.h",
"/usr/include/c++/4.7/ext/numeric_traits.h",
"/usr/include/c++/4.7/ext/type_traits.h",
"/usr/include/c++/4.7/ios",
"/usr/include/c++/4.7/iosfwd",
"/usr/include/c++/4.7/iostream",
"/usr/include/c++/4.7/istream",
"/usr/include/c++/4.7/iterator",
"/usr/include/c++/4.7/limits",
"/usr/include/c++/4.7/map",
"/usr/include/c++/4.7/new",
"/usr/include/c++/4.7/ostream",
"/usr/include/c++/4.7/set",
"/usr/include/c++/4.7/streambuf",
"/usr/include/c++/4.7/string",
"/usr/include/c++/4.7/typeinfo",
"/usr/include/c++/4.7/utility",
"/usr/include/c++/4.7/vector",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/atomic_word.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/c++allocator.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/c++config.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/c++locale.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/cpu_defines.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/ctype_base.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/ctype_inline.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/gthr-default.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/gthr.h",
"/usr/include/c++/4.7/x86_64-linux-gnu/bits/os_defines.h",
"/usr/include/ctype.h",
"/usr/include/endian.h",
"/usr/include/features.h",
"/usr/include/getopt.h",
"/usr/include/libio.h",
"/usr/include/locale.h",
"/usr/include/math.h",
"/usr/include/pthread.h",
"/usr/include/sched.h",
"/usr/include/stdio.h",
"/usr/include/stdlib.h",
"/usr/include/string.h",
"/usr/include/time.h",
"/usr/include/unistd.h",
"/usr/include/wchar.h",
"/usr/include/wctype.h",
"/usr/include/x86_64-linux-gnu/bits/byteswap.h",
"/usr/include/x86_64-linux-gnu/bits/confname.h",
"/usr/include/x86_64-linux-gnu/bits/endian.h",
"/usr/include/x86_64-linux-gnu/bits/environments.h",
"/usr/include/x86_64-linux-gnu/bits/huge_val.h",
"/usr/include/x86_64-linux-gnu/bits/huge_valf.h",
"/usr/include/x86_64-linux-gnu/bits/huge_vall.h",
"/usr/include/x86_64-linux-gnu/bits/inf.h",
"/usr/include/x86_64-linux-gnu/bits/locale.h",
"/usr/include/x86_64-linux-gnu/bits/mathcalls.h",
"/usr/include/x86_64-linux-gnu/bits/mathdef.h",
"/usr/include/x86_64-linux-gnu/bits/nan.h",
"/usr/include/x86_64-linux-gnu/bits/posix_opt.h",
"/usr/include/x86_64-linux-gnu/bits/predefs.h",
"/usr/include/x86_64-linux-gnu/bits/pthreadtypes.h",
"/usr/include/x86_64-linux-gnu/bits/sched.h",
"/usr/include/x86_64-linux-gnu/bits/select.h",
"/usr/include/x86_64-linux-gnu/bits/setjmp.h",
"/usr/include/x86_64-linux-gnu/bits/sigset.h",
"/usr/include/x86_64-linux-gnu/bits/stdio_lim.h",
"/usr/include/x86_64-linux-gnu/bits/sys_errlist.h",
"/usr/include/x86_64-linux-gnu/bits/time.h",
"/usr/include/x86_64-linux-gnu/bits/types.h",
"/usr/include/x86_64-linux-gnu/bits/typesizes.h",
"/usr/include/x86_64-linux-gnu/bits/waitflags.h",
"/usr/include/x86_64-linux-gnu/bits/waitstatus.h",
"/usr/include/x86_64-linux-gnu/bits/wchar.h",
"/usr/include/x86_64-linux-gnu/bits/wordsize.h",
"/usr/include/x86_64-linux-gnu/gnu/stubs-64.h",
"/usr/include/x86_64-linux-gnu/gnu/stubs.h",
"/usr/include/x86_64-linux-gnu/sys/cdefs.h",
"/usr/include/x86_64-linux-gnu/sys/select.h",
"/usr/include/x86_64-linux-gnu/sys/sysmacros.h",
"/usr/include/x86_64-linux-gnu/sys/types.h",
"/usr/include/xlocale.h",
"/usr/local/include/root/Buttons.h",
"/usr/local/include/root/DllImport.h",
"/usr/local/include/root/ESTLType.h",
"/usr/local/include/root/Foption.h",
"/usr/local/include/root/GuiTypes.h",
"/usr/local/include/root/Math/ParamFunctor.h",
"/usr/local/include/root/RConfig.h",
"/usr/local/include/root/RConfigure.h",
"/usr/local/include/root/RVersion.h",
"/usr/local/include/root/Riosfwd.h",
"/usr/local/include/root/Rtypeinfo.h",
"/usr/local/include/root/Rtypes.h",
"/usr/local/include/root/TArray.h",
"/usr/local/include/root/TArrayC.h",
"/usr/local/include/root/TArrayD.h",
"/usr/local/include/root/TArrayF.h",
"/usr/local/include/root/TArrayI.h",
"/usr/local/include/root/TArrayS.h",
"/usr/local/include/root/TAttAxis.h",
"/usr/local/include/root/TAttCanvas.h",
"/usr/local/include/root/TAttFill.h",
"/usr/local/include/root/TAttLine.h",
"/usr/local/include/root/TAttMarker.h",
"/usr/local/include/root/TAttPad.h",
"/usr/local/include/root/TAttText.h",
"/usr/local/include/root/TAxis.h",
"/usr/local/include/root/TBits.h",
"/usr/local/include/root/TBranch.h",
"/usr/local/include/root/TBuffer.h",
"/usr/local/include/root/TButton.h",
"/usr/local/include/root/TCanvas.h",
"/usr/local/include/root/TCanvasImp.h",
"/usr/local/include/root/TChain.h",
"/usr/local/include/root/TClass.h",
"/usr/local/include/root/TClingRuntime.h",
"/usr/local/include/root/TClonesArray.h",
"/usr/local/include/root/TCollection.h",
"/usr/local/include/root/TDataType.h",
"/usr/local/include/root/TDatime.h",
"/usr/local/include/root/TDictionary.h",
"/usr/local/include/root/TDirectory.h",
"/usr/local/include/root/TDirectoryFile.h",
"/usr/local/include/root/TError.h",
"/usr/local/include/root/TF1.h",
"/usr/local/include/root/TFile.h",
"/usr/local/include/root/TFitResultPtr.h",
"/usr/local/include/root/TFormula.h",
"/usr/local/include/root/TGenericClassInfo.h",
"/usr/local/include/root/TH1.h",
"/usr/local/include/root/TH1D.h",
"/usr/local/include/root/TH2.h",
"/usr/local/include/root/TH2D.h",
"/usr/local/include/root/THashTable.h",
"/usr/local/include/root/TIterator.h",
"/usr/local/include/root/TLatex.h",
"/usr/local/include/root/TList.h",
"/usr/local/include/root/TMap.h",
"/usr/local/include/root/TMath.h",
"/usr/local/include/root/TMathBase.h",
"/usr/local/include/root/TMatrix.h",
"/usr/local/include/root/TMatrixDBasefwd.h",
"/usr/local/include/root/TMatrixF.h",
"/usr/local/include/root/TMatrixFBasefwd.h",
"/usr/local/include/root/TMatrixFUtils.h",
"/usr/local/include/root/TMatrixFUtilsfwd.h",
"/usr/local/include/root/TMatrixFfwd.h",
"/usr/local/include/root/TMatrixT.h",
"/usr/local/include/root/TMatrixTBase.h",
"/usr/local/include/root/TMatrixTUtils.h",
"/usr/local/include/root/TMethodCall.h",
"/usr/local/include/root/TMinuit.h",
"/usr/local/include/root/TNamed.h",
"/usr/local/include/root/TObjArray.h",
"/usr/local/include/root/TObjString.h",
"/usr/local/include/root/TObject.h",
"/usr/local/include/root/TPad.h",
"/usr/local/include/root/TQObject.h",
"/usr/local/include/root/TROOT.h",
"/usr/local/include/root/TSchemaHelper.h",
"/usr/local/include/root/TSeqCollection.h",
"/usr/local/include/root/TStorage.h",
"/usr/local/include/root/TString.h",
"/usr/local/include/root/TStyle.h",
"/usr/local/include/root/TText.h",
"/usr/local/include/root/TTree.h",
"/usr/local/include/root/TUUID.h",
"/usr/local/include/root/TUrl.h",
"/usr/local/include/root/TVector2.h",
"/usr/local/include/root/TVector3.h",
"/usr/local/include/root/TVectorDfwd.h",
"/usr/local/include/root/TVectorFfwd.h",
"/usr/local/include/root/TVersionCheck.h",
"/usr/local/include/root/TVirtualPad.h",
"/usr/local/include/root/TVirtualTreePlayer.h",
"/usr/local/include/root/TVirtualX.h",
"/usr/local/include/root/snprintf.h",
"/usr/local/include/root/strlcpy.h",
"WCSimAnalysisRootDict.cc.h",
0
    };
    static const char* includePaths[] = {
"./include",
"/home/andy/work/CHIPS/WCSim/include",
"/usr/local/include/root",
"/usr/local/include/root",
"/home/andy/work/CHIPS/WCSimAnalysis/",
0
    };
    static const char* payloadCode = 
"\n"
"#ifndef ROOT_Math_VectorUtil_Cint\n"
"  #define ROOT_Math_VectorUtil_Cint\n"
"#endif\n"
"#ifndef G__VECTOR_HAS_CLASS_ITERATOR\n"
"  #define G__VECTOR_HAS_CLASS_ITERATOR\n"
"#endif\n"
"\n"
;
    static bool sInitialized = false;
    if (!sInitialized) {
      TROOT::RegisterModule("WCSimAnalysisRootDict",
        headers, allHeaders, includePaths, payloadCode,
        TriggerDictionaryInitialization_WCSimAnalysisRootDict_Impl);
      sInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_WCSimAnalysisRootDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_WCSimAnalysisRootDict() {
  TriggerDictionaryInitialization_WCSimAnalysisRootDict_Impl();
}
