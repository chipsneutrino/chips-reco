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

class WCSimEmissionProfileManager : public TObject
{
public:
	WCSimEmissionProfileManager();
	virtual ~WCSimEmissionProfileManager();
	WCSimEmissionProfiles *GetEmissionProfile(TrackType::Type type, double energy);
	WCSimEmissionProfiles *GetEmissionProfile(WCSimLikelihoodTrackBase *myTrack);
	unsigned int GetNumTracksToCache() const;
	void SetNumTracksToCache(const unsigned int numTracks);
	void CacheAtLeast(const unsigned int numTracks);

	// Get the distance a particle travels before it stops emitting photons,
	// as a smooth function of the track's energy
	double GetStoppingDistance(WCSimLikelihoodTrackBase *myTrack);

	// These wrap calls to the emission profile object:
	TH1F *GetEnergyHist(const TrackType::Type &type);
	bool EnergiesInSameBin(const TrackType::Type &type, const double energy1, const double energy2);
	TH1F *GetRho(const TrackType::Type &type, const double &energy);
	std::pair<TH2F *, TH2F *> GetG(const TrackType::Type &type, const double &energy);
	double GetLightFlux(WCSimLikelihoodTrackBase *myTrack);
	double GetTrackLengthForPercentile(WCSimLikelihoodTrackBase *myTrack, double percentile);
	double GetTrackLength(const TrackType::Type &type, const double &energy);
	std::vector<Double_t> GetRhoIntegrals(std::vector<Int_t> sPowers, WCSimLikelihoodTrackBase *myTrack,
										  Double_t startS, Double_t endS, Bool_t multiplyByWidth = kTRUE);
	std::vector<double> GetRhoGIntegrals(WCSimLikelihoodTrackBase *myTrack, WCSimLikelihoodDigit *myDigit,
										 std::vector<int> sPowers, double cutoffS, bool multiplyByWidth = true);

	std::vector<double> GetFourNearestEnergies(WCSimLikelihoodTrackBase *myTrack);
	std::vector<double> GetFourNearestEnergies(const TrackType::Type &type, const double &energy);
	std::vector<TH1F *> GetFourNearestSForTimeHists(WCSimLikelihoodTrackBase *myTrack);
	std::vector<TH2F *> GetFourNearestSCosThetaForTimeHists(WCSimLikelihoodTrackBase *myTrack);

private:
	void AddNewProfile(const TrackType::Type &type, const double &energy);
	std::vector<double> fEnergies;

	void ResetSForTimeHists();
	void ResetSCosThetaForTimeHists();
	std::vector<TH1F *> fSForTimeHists;
	std::vector<TH2F *> fSCosThetaForTimeHists;
	int fLastNearbyEnergyBin;
	TrackType::Type fLastNearbyType;

	unsigned int fNumTracksToCache;
	std::map<TrackType::Type, std::map<double, WCSimEmissionProfiles *>> fProfileMap;
	std::map<TrackType::Type, std::vector<double>> fOrderOfUse;
	void UpdateRecentlyUsed(const TrackType::Type &type, const double &energy);
	double GetLeastRecentlyUsedEnergy(const TrackType::Type &type);
	double GetMostRecentlyUsedEnergy(const TrackType::Type &type);
	void RemoveFromRecentEnergies(const TrackType::Type &type, const double &energy);
	unsigned int GetNumCached(const TrackType::Type &type);

	TrackType::Type fLastLengthType;
	double fLastLengthEnergy;
	double fLastLength;

	ClassDef(WCSimEmissionProfileManager, 1)
};

#endif
