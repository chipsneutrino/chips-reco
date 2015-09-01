/*
 * WCSimPiZeroFitter.cc
 *
 *  Created on: 21 Aug 2015
 *      Author: ajperch
 */
#include "WCSimAnalysisConfig.hh"
#include "WCSimChargePredictor.hh"
#include "WCSimFitterConfig.hh"
#include "WCSimFitterParameters.hh"
#include "WCSimFitterPlots.hh"
#include "WCSimFitterTree.hh"
#include "WCSimFitterInterface.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodDigitArray.hh"
#include "WCSimLikelihoodFitter.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimReco.hh"
#include "WCSimRecoSeed.hh"
#include "WCSimRecoDigit.hh"
#include "WCSimRecoEvent.hh"
#include "WCSimRecoFactory.hh"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimTimeLikelihood.hh"
#include "WCSimTotalLikelihood.hh"
#include "WCSimTrackParameterEnums.hh"

#include "TClonesArray.h"
#include "TCollection.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TMinuit.h"
#include "TPrincipal.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include <string>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include <map>
#include <algorithm>

#include "WCSimPiZeroFitter.hh"

#include "TCanvas.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVectorD.h"


#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPiZeroFitter)
#endif

bool RingSort(const std::pair<WCSimRecoRing*,double> &a, const std::pair<WCSimRecoRing*,double> &b);

/**
 * @todo Remove hardcoding of track type
 */
WCSimPiZeroFitter::WCSimPiZeroFitter()
{
  fFitterPlots = NULL;
  fFitterTree = 0x0;
  fTotalLikelihood = NULL;
  fRootEvent = NULL;
  fLikelihoodDigitArray = NULL;
  fTrueLikelihoodTracks = NULL;
  fCalls = 0;

  fUseHoughFitterForSeed = WCSimAnalysisConfig::Instance()->GetUseHoughFitterForSeed();
  ResetEvent();
}

WCSimPiZeroFitter::~WCSimPiZeroFitter()
{
}

void WCSimPiZeroFitter::FitAfterFixingDirectionAndEnergy(
		std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> pair)
{
	SetStartingTracks(pair);
	FreeVertex();
	FreeEnergy();
	FreeDirection();

	FixDirection();
	FixEnergy();
	Fit("Simplex");

	FreeDirection();
	FreeEnergy();
	Fit("Simplex");
}

void WCSimPiZeroFitter::FitAfterFixingEnergy(
		std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> pair)
{
	SetStartingTracks(pair);
	FreeVertex();
	FreeDirection();

	FixEnergy();
	Fit("Simplex");
	FreeEnergy();

	Fit("Simplex");
}

void WCSimPiZeroFitter::SetStartingTracks(
		std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> tracks)
{

    // Set first track starting values
    if(CanSetParam(0, FitterParameterType::kVtxX)){
        fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxX, tracks.first->GetX());
    }
    if(CanSetParam(0, FitterParameterType::kVtxY)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxY, tracks.first->GetY());
    }
    if(CanSetParam(0, FitterParameterType::kVtxZ)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxZ, tracks.first->GetZ());
    }
    if(CanSetParam(0, FitterParameterType::kVtxT)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxT, tracks.first->GetT());
    }
    if(CanSetParam(0, FitterParameterType::kDirTh)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kDirTh, tracks.first->GetTheta());
    }
    if(CanSetParam(0, FitterParameterType::kDirPhi)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kDirPhi, tracks.first->GetPhi());
    }
    if(CanSetParam(0, FitterParameterType::kEnergy)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kEnergy, tracks.first->GetE());
    }
    if(CanSetParam(0, FitterParameterType::kConversionDistance)){
            fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kConversionDistance, tracks.first->GetConversionDistance());
    }

    // Set second track starting values
    if(CanSetParam(1, FitterParameterType::kVtxX)){
        fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxX, tracks.second->GetX());
    }
    if(CanSetParam(1, FitterParameterType::kVtxY)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxY, tracks.second->GetY());
    }
    if(CanSetParam(1, FitterParameterType::kVtxZ)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxZ, tracks.second->GetZ());
    }
    if(CanSetParam(1, FitterParameterType::kVtxT)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kVtxT, tracks.second->GetT());
    }
    if(CanSetParam(1, FitterParameterType::kDirTh)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirTh, tracks.second->GetTheta());
    }
    if(CanSetParam(1, FitterParameterType::kDirPhi)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirPhi, tracks.second->GetPhi());
    }
    if(CanSetParam(1, FitterParameterType::kEnergy)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kEnergy, tracks.second->GetE());
    }
    if(CanSetParam(1, FitterParameterType::kConversionDistance)){
            fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kConversionDistance, tracks.second->GetConversionDistance());
    }

    return;

}

bool WCSimPiZeroFitter::CanSetParam(const unsigned int &iTrack, const FitterParameterType::Type &type)
{
    // Only see the parameters if they are not requested to be fixed:
    // Also check if any of the parameters are joined with any other tracks: we should
    // let the lower-numbered track win in this case, because it has the higher Hough peak
	return ((!fFitterTrackParMap.GetIsFixed(iTrack, type))
			 && (WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(iTrack, type) >= iTrack));

}

UInt_t WCSimPiZeroFitter::GetNPars()
{
  // Do we know how to fit this number of tracks?
  unsigned int nPars = WCSimFitterConfig::Instance()->GetNumIndependentParameters();
  return nPars;
}

void WCSimPiZeroFitter::RunFits()
{
	std::cout << "WCSimPiZeroFitter::RunFits()" << std::endl;
	FitEventNumber(0);

}

void WCSimPiZeroFitter::FitEventNumber(Int_t iEvent) {
	std::cout << "WCSimPiZeroFitter::FitEventNumber " << iEvent << std::endl;
	fIsFirstCall = true;
	fEvent = iEvent;
	fCalls = 0;
	WCSimInterface::Instance()->BuildEvent(iEvent);
	fRootEvent = WCSimInterface::Instance()->GetWCSimEvent(iEvent);
	fLikelihoodDigitArray = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray(iEvent);
	std::cout << "There are " << fLikelihoodDigitArray->GetNDigits() << " digits" << std::endl;
	fFitterTrackParMap.Set();


	if(fTotalLikelihood != NULL)
	{
		delete fTotalLikelihood;
		fTotalLikelihood = NULL;
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}
	else
	{
		fTotalLikelihood = new WCSimTotalLikelihood(fLikelihoodDigitArray);
	}

	//////////////////////////////////////////////////////////
	// Build all the combinations of tracks to try as seeds:
	//////////////////////////////////////////////////////////

	// First we fit a single electron track using the Hough transform.  This gets us a vertex
	// and a direction to use as the jumping off point to generate the seeds for the fitter
	SeedWithSingleElectron();

	// Try various combinations of conversion distance for each of the two photons.  We'll run the fitter
	// several times with different track directions for each of these combinations
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > allSeedTracks;
	allSeedTracks = IterateOverConversionDistances();
		/* |-> Fixes the conversion distance of each track to a pre-chosen value, and tries all combinations
		 	|-> Each combination iterates over the direction of the single track seed by perturbing the direction slightly
		     |--> For each perturbation of the single track, does a grid search to find the best
		 	 	  directions for the second track and saves one where the tracks are close in energy and
				  one where they're separate
		  Returns a vector of all the track pairs we want to try for the seed
		*/


	////////////////////////////////////////////////////////////////////////////////////
	// Now try all the seed tracks one at a time in the fitter, and pick the best result
	////////////////////////////////////////////////////////////////////////////////////

	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> >::iterator seedTrackItr;
	seedTrackItr = allSeedTracks.begin();
	std::cout << "There are " << allSeedTracks.size() << " seeds" << std::endl;

	// We try two minimiser algorithms:
	// 1. Fix the energy and direction of both tracks, run the fit to find a minimum
	//    Then free everything and re-run the fit, starting at this minimum
	// 2. As above, except it only fixes the energies and leaves the directions free
	while(seedTrackItr != allSeedTracks.end())
	{
		// These will update fMinimum and the fBestFits if the minimiser improves on the previous best
		FitAfterFixingDirectionAndEnergy(*seedTrackItr);
		FitAfterFixingEnergy(*seedTrackItr);
		++seedTrackItr;
	}

	fTrueLikelihoodTracks = WCSimInterface::Instance()->GetTrueLikelihoodTracks();
    if( fMinimum > 0 )
    {
		  FillPlots();
		  FillTree();
		  FillHitComparison();
    }
    else
    {
      fFitterTree->FillRecoFailures(iEvent);
    }

  	std::cout << "Fitted event number " << iEvent << std::endl;
	return;
}

std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > WCSimPiZeroFitter::IterateOverConversionDistances()
{
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > allTracks;

	static const int numDists = 2;
	static const double firstDist = 50.0; //cm
	static const double secondDist = 250.0; //cm
	double conversionDistancesToTry[2] = { firstDist, secondDist };
	for(int iTrack1Dist = 0; iTrack1Dist < numDists; ++iTrack1Dist )
	{
		for(int iTrack2Dist = 0; iTrack2Dist < numDists; ++iTrack2Dist)
		{
			std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > trackOptions;
			trackOptions = IterateOverFirstTrackPerturbations(conversionDistancesToTry[iTrack1Dist],
															  conversionDistancesToTry[iTrack2Dist]);
			for(unsigned int iTrack = 0; iTrack < trackOptions.size(); ++iTrack)
			{
				allTracks.push_back(trackOptions.at(iTrack));
			}
		}
	}

	return allTracks;
}

std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > WCSimPiZeroFitter::IterateOverFirstTrackPerturbations(const double &convDistTrack1, const double &convDistTrack2)
{
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > allTracks;

	std::vector<TVector3> directionsToTry = GetFirstTrackDirectionsToTry();
	std::vector<TVector3>::iterator directionItr = directionsToTry.begin();
	while(directionItr != directionsToTry.end())
	{
		std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase *> > tracks;

		std::cout << "I think the distance is " << std::distance(directionsToTry.begin(), directionItr) <<std::endl;

		tracks = GridSearchOverSecondTrackDirection(
				convDistTrack1, convDistTrack2,
				*directionItr,
				std::distance(directionsToTry.begin(), directionItr));

		// Append the new track vector to allTracks, after checking that tracks
		// isn't empty (it could be if there's a big energy imbalance between the
		// two tracks  given this separation angle)
		if(tracks.size() > 0)
		{
			allTracks.push_back(tracks.at(0));
		}
		std::cout << "Size of vector returned by grid search = " << tracks.size() << std::endl;
		++directionItr;
	}

	return allTracks;
}

std::vector<TVector3> WCSimPiZeroFitter::GetFirstTrackDirectionsToTry()
{
	// Use the single track seed to define the plane containing its vertex, normal to its direction:
	TVector3 planeContains = fSingleElectronSeed.fVertex;
	TVector3 planeUnitNormal = (fSingleElectronSeed.fDirection).Unit();
	std::cout << "Point in plane = " << std::endl;
	planeContains.Print();
	std::cout << "Normal to plane = " << std::endl;
	planeUnitNormal.Print();


	TVector3 planeAxis1, planeAxis2; // Orthogonal unit vectors in the plane
	planeAxis1 = planeUnitNormal.Cross(TVector3(0,1,0));
	if(planeAxis1.Mag() < 1e-10)
	{
		// If our plane normal happens to be exactly along the y-axis
		// then the above cross product will be zero and we should use
		// a different axis to generate the first vector in the plane
		planeAxis1 = planeUnitNormal.Cross(TVector3(1,0,0));
	}
	planeAxis1 = planeAxis1.Unit();
	planeAxis2 = planeAxis1.Cross(planeUnitNormal).Unit();

	WCSimLikelihoodDigitArray * myLDA = WCSimInterface::Instance()->GetWCSimLikelihoodDigitArray();
	int nDigits = myLDA->GetNDigits();

	std::vector<TVector2> projections;
	TCanvas * can = new TCanvas("can","",800,600);
	TGraph * gr = new TGraph();
	gr->SetName("gr");

	// We can't just take the component of the PMT position relative to the track vertex that
	// is perpendicular to the track direction because our ring might straddle an edge, e.g. between
	// the barrel and endcap.  For example, with a track directed along x, everything on the roof
	// would have the same projection
	// Instead we need to project all the PMTs into a plane the same distance away from the vertex
	// in the direction parallel to the track momentum.  Because we are ultimately only interested in the
	// direction, this can be any distance so long as it's constant - I've picked the cylinder radius
	// because it's a convenient scale to visualise
	double projectToDistance = WCSimGeometry::Instance()->GetCylRadius();

	for(int iDigit = 0; iDigit < nDigits; ++iDigit)
	{
		if(myLDA->GetDigit(iDigit)->GetQ() < 1.0) { continue; }
		TVector3 digitPos = myLDA->GetDigit(iDigit)->GetPos();
		double distanceInTrackDir = (digitPos - planeContains).Dot(planeUnitNormal);
		// Move the PMT out so they all sit in the same plane
		if(fabs(projectToDistance/distanceInTrackDir) > 10) { continue; }
		TVector3 extendedDigitPos = digitPos * fabs(projectToDistance/distanceInTrackDir);

		// Then project back into the plane containing the vertex
		TVector3 projection = extendedDigitPos - ((extendedDigitPos - planeContains).Dot(planeUnitNormal)) * planeUnitNormal;

		// Project the points into an (x',y') plane, where x' = planeAxis1 and y' = planeAxis2
		projections.push_back(TVector2(projection.Dot(planeAxis1), projection.Dot(planeAxis2)));
	}

	// Now calculate the covariance matrix of these points in the (x',y') plane
	// Turns out there's a ROOT class that can do this for us
	TPrincipal myPrincipal(2, ""); // Args: nVariables, (N)ormalise cov matrix and store (D)ata
	double meanX = 0.0;
	double meanY = 0.0;
	for(unsigned int iProj = 0; iProj < projections.size(); ++iProj)
	{
		double point[2] = {projections.at(iProj).X(), projections.at(iProj).Y() };
		gr->SetPoint(gr->GetN(), point[0], point[1]);
		myPrincipal.AddRow(point);
		meanX += point[0]/projections.size();
		meanY += point[1]/projections.size();
	}

	gr->Draw("AP");
	can->SaveAs("can1.C");

	// Calculate the eigenvectors of the covariance matrix
	// to get the major and minor axes of the covariance ellipse

	myPrincipal.GetCovarianceMatrix()->Print(); // This only stores the lower half of the matrix by default
	double scaleEigenvalue = (*(myPrincipal.GetCovarianceMatrix()))(0,0);
	myPrincipal.MakePrincipals();  // This step performs a normalisation of the covariance matrix and calculates the upper half
	scaleEigenvalue /= (*(myPrincipal.GetCovarianceMatrix()))(0,0);  // Multiply the eigenvalue by this factor to undo the scaling
										  	  	  	  	  	  	  	 // All this does is make the plot look nicer; it's not vital

	const TVectorD eigenvalues = *(myPrincipal.GetEigenValues());
	const TMatrixD eigenvectors = *(myPrincipal.GetEigenVectors());

	// This will be 1 if the second eigenvector has a larger eigenvalue than the first
	// So it can be used as the index for the eigenvector corresponding to the major axis
	bool secondVectorIsMajorAxis = (eigenvalues[1] > eigenvalues[0]);
	// Eigenvectors should be arranged in eigenvalue order, but this makes sure

	// Work out the major and minor axes of the ellipse in the normal (x,y,z) coordinates of the detector
	TVector2 eigenvector(eigenvectors(secondVectorIsMajorAxis,0), eigenvectors(secondVectorIsMajorAxis,1));
	double phi = eigenvector.Phi();  // Angle to rotate th elipse by - again, just for the plot
	TVector3 majorAxis = eigenvector.X() * planeAxis1 + eigenvector.Y() * planeAxis2;
	eigenvector = TVector2(eigenvectors(!secondVectorIsMajorAxis,0), eigenvectors(!secondVectorIsMajorAxis,1));
	TVector3 minorAxis = eigenvector.X() * planeAxis1 + eigenvector.Y() * planeAxis2;

	// Normalise them
	majorAxis = majorAxis.Unit();
	minorAxis = minorAxis.Unit();

	// Draw an ellipse
	TEllipse el(meanX, meanY,
				TMath::Sqrt(2.30*(eigenvalues[0]) * scaleEigenvalue),
				TMath::Sqrt(2.30*(eigenvalues[1]) * scaleEigenvalue),
				0,360,
				phi);
	el.SetLineColor(kBlue);
	el.SetFillStyle(0);
	el.SetLineWidth(2);
	el.Draw("");
	can->SaveAs("withEllipse.C");
	can->SaveAs("withEllipse.png");



	// Now work out the other track (theta, phi) pairings to try by rotating the seed track about the
	// semimajor and semiminor axes
	std::vector<TVector3> directionsToTry;  // The final set of unit vectors corresponding to track 1 directions

	static const int numRotMaj = 7; // Rotations about the major axis
	double rotationsAboutMajor[numRotMaj] = {  0.0,
											   0.05 * M_PI, -0.05 * M_PI,
											   0.45 	     , -0.45,
											   0.20 * M_PI, -0.20 * M_PI };
	for(int i = 0; i < numRotMaj; ++i)
	{
		TVector3 tmp = fSingleElectronSeed.fVertex; // TVector3::Rotate overwrites the vector, so need a copy
		tmp.Rotate(rotationsAboutMajor[i], majorAxis);
		directionsToTry.push_back(tmp);
	}

	static const int numRotMin = 2; // Rotations about the minor axis
	double rotationsAboutMinor[numRotMin] = { 0.05*M_PI, -0.05 * M_PI };
	for(int jRot = 0; jRot < numRotMin; ++jRot)
	{
		TVector3 tmp = fSingleElectronSeed.fVertex;
		tmp.Rotate(rotationsAboutMinor[jRot], minorAxis);
		directionsToTry.push_back(tmp);
	}
	delete can;
	delete gr;
	return directionsToTry;
}

std::vector<std::pair<WCSimLikelihoodTrackBase *, WCSimLikelihoodTrackBase *> > WCSimPiZeroFitter::GridSearchOverSecondTrackDirection(
		const double &convDistTrack1,
		const double &convDistTrack2,
		const TVector3 &track1Dir,
		int numTimesAlready)
{
	// std::cout << "*** WCSimPiZeroFitter::GridSearchOverSecondTrackDirection *** " << std::endl;

	const unsigned int nBins = 5; // Number of bins in theta and phi
	assert(nBins);

	// Make first track:
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxX, fSingleElectronSeed.fVertex.X());
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxY, fSingleElectronSeed.fVertex.Y());
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxZ, fSingleElectronSeed.fVertex.Z());
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kVtxT, fSingleElectronSeed.fTime);
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kEnergy, fSingleElectronSeed.fEnergy);
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kConversionDistance, convDistTrack1);

	double track1Theta = track1Dir.Theta();
	double track1Phi = track1Dir.Phi();
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kDirTh, track1Theta);
	fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kDirPhi, track1Phi);

	if(   WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(1, FitterParameterType::kVtxX) != 0
	   || WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(1, FitterParameterType::kVtxY) != 0
	   || WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(1, FitterParameterType::kVtxZ) != 0
	   || WCSimFitterConfig::Instance()->GetTrackIsJoinedWith(1, FitterParameterType::kVtxT) != 0)
	{
		std::cerr << "Error: trying to fit a pi0 by the two electrons don't have the same vertex" << std::endl;
		assert(0);
	}

	// Default the best values to something junk so we can tell when they've been set
	double best2LnL = -999.9;
	double bestTheta = -999.9;
	double bestPhi = -999.9;
	double bestEnergy = -999.9;

	// Start the grid search
	fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kConversionDistance, convDistTrack2);
	for(unsigned int thetaBin = 0; thetaBin < nBins; ++thetaBin)
	{
		double theta = thetaBin * M_PI/nBins;
		fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirTh, theta);

		for(unsigned int phiBin = 0; phiBin < nBins; ++phiBin)
		{
			double phi = phiBin * 2.0 * M_PI/nBins;
			std::cout << "Grid search number " << numTimesAlready+1 << "  Bin " << phiBin+1 + thetaBin * nBins << "/" << nBins*nBins << "   (theta, phi) = (" << theta << ", " << phi << ")" << std::endl;
			fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirPhi, phi);

			TVector3 track2Dir;
			track2Dir.SetMagThetaPhi(1, theta, phi);
			double energy = GetPiZeroSecondTrackEnergy(track1Dir, fSingleElectronSeed.fEnergy, track2Dir);
			if(energy > 5000)
			{
				continue;
			}
			if(   numTimesAlready > 0
			   && (   energy < 0.05 * fSingleElectronSeed.fEnergy
				   || fSingleElectronSeed.fEnergy < 0.05*energy )
			){ continue; }
			fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kEnergy, energy);

			std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
			Double_t * x = &(startVals[0]);
     		double tempMin = WrapFunc(x);

     		if( tempMin < best2LnL || best2LnL < 0 )
			{
				best2LnL = tempMin;
				bestTheta = theta;
				bestPhi = phi;
				bestEnergy = energy;
			}
		}
	}

	// Return the best-fit
	std::vector<std::pair<WCSimLikelihoodTrackBase*, WCSimLikelihoodTrackBase*> > best;
	if( best2LnL > 0)
	{
		std::map<FitterParameterType::Type, double> extraPars;
		extraPars[FitterParameterType::kConversionDistance] = convDistTrack1;
		WCSimLikelihoodTrackBase * firstTrack = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::PhotonLike,
				fSingleElectronSeed.fVertex.X(), fSingleElectronSeed.fVertex.Y(), fSingleElectronSeed.fVertex.Z(),
				fSingleElectronSeed.fTime,
				track1Theta, track1Phi,
				fSingleElectronSeed.fEnergy,
				extraPars);


		std::map<FitterParameterType::Type, double> extraPars2;
		extraPars[FitterParameterType::kConversionDistance] = convDistTrack2;
		WCSimLikelihoodTrackBase * secondTrack = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::PhotonLike,
				fSingleElectronSeed.fVertex.X(), fSingleElectronSeed.fVertex.Y(), fSingleElectronSeed.fVertex.Z(),
				fSingleElectronSeed.fTime,
				bestTheta,
				bestPhi,
				bestEnergy,
				extraPars2);
		best.push_back(std::make_pair(firstTrack, secondTrack));
	}
	std::cout << "Best 2LnL was " << best2LnL << std::endl;
	std::cout << "Size of best-fit vector was " << best.size() << std::endl;
	return best;
}


// Get the energy of photon track 1 given the energy of photon track 0
// assuming they arise from the decay of a pi zero (and so have its invariant
// mass)
// \param x The array of current values that will be passed to/returned from WrapFuncPiZero
//          (which comes from the FitterTrackParMap)
Double_t WCSimPiZeroFitter::GetPiZeroSecondTrackEnergy(const TVector3 &track1Dir, const double &track1Energy, const TVector3 &track2Dir )
{
  double cosThTrk1 = track1Dir.CosTheta();
  double sinThTrk1 = sin(track1Dir.Theta());
  double cosPhiTrk1 = cos(track1Dir.Phi());
  double sinPhiTrk1 = sin(track1Dir.Phi());
  double cosThTrk2 = track2Dir.CosTheta();
  double sinThTrk2 = sin(track2Dir.Theta());
  double cosPhiTrk2 = cos(track2Dir.Phi());
  double sinPhiTrk2 = sin(track2Dir.Phi());


  double massOfPiZeroMeV = 134.9766; // PDG 2014: http://pdg.lbl.gov/2014/listings/rpp2014-list-pi-zero.pdf

  double cosTheta01 =   sinThTrk1 * cosPhiTrk1 * sinThTrk2 * cosPhiTrk2
                      + sinThTrk1 * sinPhiTrk1 * sinThTrk2 * sinPhiTrk2
                      + cosThTrk1 * cosThTrk2;
                      // This is the dot product of unit vectors in the directions of track 0 and track 1
                      // i.e. the cosine of the angle between the two photons
  double track2Energy = 0.5 * massOfPiZeroMeV * massOfPiZeroMeV / (track1Energy * (1.0 - cosTheta01));

  TLorentzVector fourMom0( track1Energy * sinThTrk1 * cosPhiTrk1, track1Energy * sinThTrk1 * sinPhiTrk1, track1Energy * cosThTrk1, track1Energy);
  TLorentzVector fourMom1( track2Energy * sinThTrk2 * cosPhiTrk2, track2Energy * sinThTrk2 * sinPhiTrk2, track2Energy * cosThTrk2, track2Energy);
  // std::cout << "E1 = " << track1Energy << "   E2 = " << track2Energy << "   Separation = " << TMath::ACos(cosTheta01) << "  Total four-momentum = " << (fourMom0 + fourMom1).Mag() << std::endl;

  return track2Energy;

}

void WCSimPiZeroFitter::SeedWithSingleElectron()
{

	std::cout << " *** WCSimPiZeroFitter::SeedWithSingleElectron() *** " << std::endl;

	// Run the old Hough transform reco
	WCSimRecoSeed * myReco = dynamic_cast<WCSimRecoSeed*>(WCSimRecoFactory::Instance()->MakeReco("seed")); // This calls new() - delete when done
	myReco->SetNumberOfTracks(1);

	WCSimRecoEvent* recoEvent = WCSimInterface::RecoEvent();
	std::vector<WCSimRecoEvent*> slicedEvents = myReco->RunSeed(recoEvent);

	// Make a vector of all of the available rings we have.
	std::vector<WCSimRecoRing*> ringVec;
	std::vector<std::pair<WCSimRecoRing*, double> > otherRings;
	std::vector<double> ringTime;
	// Add the primary ring from each slice first.
	for (unsigned int e = 0; e < slicedEvents.size(); ++e)
	{
		//ringVec.push_back(slicedEvents[e]->GetRing(0));
		//ringTime.push_back(slicedEvents[e]->GetVtxTime());
		for (int r = 0; r < slicedEvents[e]->GetNRings(); ++r)
		{
			std::pair<WCSimRecoRing*, double> tempPair;
			tempPair = std::make_pair<WCSimRecoRing*, double>(
					slicedEvents[e]->GetRing(r), slicedEvents[e]->GetVtxTime());
			otherRings.push_back(tempPair);
		}
	}
	// Sort the otherRings vector by height and add to the main vector
	std::sort(otherRings.begin(), otherRings.end(), RingSort);
	for (unsigned int r = 0; r < otherRings.size(); ++r)
	{
		ringVec.push_back(otherRings[r].first);
		ringTime.push_back(otherRings[r].second);
	}

	// Now we need to loop over all the track parameters available to the fitter
	// and set them to the corresponding seed parameter
	double seedX, seedY, seedZ, seedT;
	double dirX, dirY, dirZ;
	seedX = ringVec[0]->GetVtxX();
	seedY = ringVec[0]->GetVtxY();
	seedZ = ringVec[0]->GetVtxZ();
	seedT = ringTime[0];
	dirX = ringVec[0]->GetDirX();
	dirY = ringVec[0]->GetDirY();
	dirZ = ringVec[0]->GetDirZ();

	fSingleElectronSeed.fVertex = TVector3(seedX, seedY, seedZ);
	fSingleElectronSeed.fDirection = TVector3(dirX, dirY, dirZ);
	fSingleElectronSeed.fTime = seedT;
	FitSingleElectronSeedEnergy();
	delete myReco;

	// Need to delete the elements of slicedEvents as they are not destroyed by WCSimRecoSlicer
	for (unsigned int v = 0; v < slicedEvents.size(); ++v)
	{
		delete slicedEvents[v];
	}
}

void WCSimPiZeroFitter::FitSingleElectronSeedEnergy()
{
	std::cout << "*** WCSimPiZeroFitter::FitSingleElectronSeedEnergy *** " << std::endl;
	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

	// Alternatively: use a different algorithm to check the minimizer works
	// ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "GSLSimAn");


	min->SetMaxFunctionCalls(500);
	min->SetMaxIterations(1);
	min->SetPrintLevel(3);
	//min->SetTolerance(0.);
	min->SetErrorDef(1.0);
	min->SetStrategy(2);
	std::cout << " Tolerance = " << min->Tolerance() << std::endl;

	// Tell the minimizer the function to minimize
	// We have to wrap it in this functor because it's a member function of this class
	const unsigned int nPars = 1;
	ROOT::Math::Functor func(this,&WCSimPiZeroFitter::WrapFuncSeedTrackEnergy, nPars);


	// Tell the minimizer the functor and variables to consider
	min->SetFunction(func);
	double minVal = fFitterTrackParMap.GetMinValue(0, FitterParameterType::kEnergy);
	double maxVal = fFitterTrackParMap.GetMaxValue(0, FitterParameterType::kEnergy);
	double step = fFitterTrackParMap.GetStep(0, FitterParameterType::kEnergy);
	min->SetLimitedVariable(0, "Energy", 0, step, minVal, maxVal);
	min->Minimize();

	double bestEnergy = min->X()[0];
	fSingleElectronSeed.fEnergy = GetFirstPhotonEnergyEstimate(bestEnergy);
	std::cout << "Best-fit single electron energy = " << bestEnergy
			  << " -> first photon energy estimate = " << fSingleElectronSeed.fEnergy << std::endl;
	return;

}

double WCSimPiZeroFitter::WrapFuncSeedTrackEnergy(const double * x)
{
	std::vector<WCSimLikelihoodTrackBase *> tracksToFit;
	WCSimLikelihoodTrackBase * track = WCSimLikelihoodTrackFactory::MakeTrack(
			TrackType::ElectronLike,
			fSingleElectronSeed.fVertex.X(),
			fSingleElectronSeed.fVertex.Y(),
			fSingleElectronSeed.fVertex.Z(),
			fSingleElectronSeed.fTime,
			fSingleElectronSeed.fDirection.Theta(),
			fSingleElectronSeed.fDirection.Phi(),
			x[0]);
	tracksToFit.push_back(track);

  Double_t minus2LnL = 0.0;
  fTotalLikelihood->SetTracks(tracksToFit);
  minus2LnL = fTotalLikelihood->Calc2LnL();

  std::cout << "Function evaluated " << ++fCalls << " times for event " << fEvent << " ... " << " -2Ln(L) = " << minus2LnL << std::endl << std::endl;
  return minus2LnL;
}


double WCSimPiZeroFitter::GetFirstPhotonEnergyEstimate(
		const double& singleElectronEnergy)
{
	return singleElectronEnergy;
}
