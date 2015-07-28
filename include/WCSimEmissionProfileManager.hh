#ifndef WCSIMEMISSIONPROFILEMANAGER_HH
#define WCSIMEMISSIONPROFILEMANAGER_HH

#include "WCSimTrackParameterEnums.hh"
#include "WCSimEmissionProfiles.hh"
#include <TObject.h>
#include <map>
#include <vector>


class TH1F;
class TH2F;
class WCSimLikelihoodTrackBase;

class WCSimEmissionProfileManager : public TObject {
  public:
      WCSimEmissionProfileManager();
      virtual ~WCSimEmissionProfileManager();
      WCSimEmissionProfiles * GetEmissionProfile(TrackType::Type type, double energy);
      WCSimEmissionProfiles * GetEmissionProfile(WCSimLikelihoodTrackBase * myTrack);
      unsigned int GetNumTracksToCache() const;
      void SetNumTracksToCache(const unsigned int numTracks);
      void CacheAtLeast(const unsigned int numTracks);

      // These wrap calls to the emission profile object:
      double GetStoppingDistance(WCSimLikelihoodTrackBase * myTrack);
      TH1F * GetEnergyHist(const TrackType::Type &type);
      TH1F * GetRho(const TrackType::Type &type, const double &energy);
      std::pair<TH2F *, TH2F *>  GetG(const TrackType::Type &type, const double &energy);
      double GetLightFlux(WCSimLikelihoodTrackBase * myTrack);
      double GetTrackLengthForPercentile(WCSimLikelihoodTrackBase * myTrack, double percentile);
	    std::vector<Double_t> GetRhoIntegrals(std::vector<Int_t> sPowers, WCSimLikelihoodTrackBase* myTrack, Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
      std::vector<double> GetRhoGIntegrals(WCSimLikelihoodTrackBase* myTrack, WCSimLikelihoodDigit* myDigit, std::vector<int> sPowers, double cutoffS, bool multiplyByWidth = true);


  private:
      void AddNewProfile(const TrackType::Type &type, const double &energy);

      unsigned int fNumTracksToCache;
      std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles*> > fProfileMap;
      std::map<TrackType::Type, std::vector<double> > fOrderOfUse;
      void UpdateRecentlyUsed(const TrackType::Type &type, const double &energy);
      double GetLeastRecentlyUsedEnergy(const TrackType::Type &type);
      double GetMostRecentlyUsedEnergy(const TrackType::Type &type);
      void RemoveFromRecentEnergies(const TrackType::Type &type, const double &energy);
      unsigned int GetNumCached(const TrackType::Type &type);
  
  ClassDef(WCSimEmissionProfileManager,1)
};

#endif
