#include "WCSimDisplayFactory.hh"
#include "WCSimDisplay.hh"
#include "WCSimDisplayAB.hh"
#include "WCSimEveDisplay.hh"
#include "WCSimVertexViewer.hh"

#include <string>
#include <iostream>
#include <cassert>

ClassImp(WCSimDisplayFactory)

static WCSimDisplayFactory* fgDisplayFactory = 0;

WCSimDisplayFactory* WCSimDisplayFactory::Instance()
{
  if(!fgDisplayFactory){
    fgDisplayFactory = new WCSimDisplayFactory();
  }

  // die if finder hasn't actually been created
  if(!fgDisplayFactory){
    assert(fgDisplayFactory);
  }

  // can do re-setting here
  if(fgDisplayFactory){

  }

  return fgDisplayFactory;
}

WCSimDisplay* WCSimDisplayFactory::MakeDisplay(const char *name)
{
  std::cout << " *** WCSimDisplayFactory::MakeDisplay(...) *** " << std::endl;

  if( strcmp(name,"AB")==0 
   || strcmp(name,"default")==0 ) {
    if( fDisplayAB==0 ){
      std::cout << "   Building Display: " << name << std::endl;
      fDisplayAB = new WCSimDisplayAB();
    }
    std::cout << "   Returning Display: " << name << std::endl;
    return fDisplayAB;
  }  

  else if( strcmp(name,"Eve")==0 
        || strcmp(name,"EVE")==0 
        || strcmp(name,"3D")==0 ) {
    if( fDisplayEve==0 ){
      std::cout << "   Building Display: " << name << std::endl;
      fDisplayEve = new WCSimEveDisplay();
    }
    std::cout << "   Returning Display: " << name << std::endl;
    return fDisplayEve;
  }

  else if( strcmp(name,"Vtx")==0 
        || strcmp(name,"VTX")==0
        || strcmp(name,"Vertex")==0 ) { 
    if( fDisplayVtx==0 ){
      std::cout << "   Building Display: " << name << std::endl;
      fDisplayVtx = new WCSimVertexViewer();
    }
    std::cout << "   Returning Display: " << name << std::endl;
    return fDisplayVtx;
  }

  else {
    std::cout << "   <warning> unknown display: " << name << std::endl;
  }

  return 0;
}

WCSimDisplayFactory::WCSimDisplayFactory()
{
  // list of possible displays
  fDisplayAB = 0;
  fDisplayEve = 0;
  fDisplayVtx = 0;
}

WCSimDisplayFactory::~WCSimDisplayFactory()
{
  if( fDisplayAB )  delete fDisplayAB;
  if( fDisplayEve ) delete fDisplayEve;
  if( fDisplayVtx ) delete fDisplayVtx;
}
