/*
 * WCSimPiZeroElectronAdjuster.cc
 *
 *  Created on: 9 Sep 2015
 *      Author: ajperch
 */

#include "WCSimPiZeroElectronAdjuster.hh"

#include "WCSimFitterConfig.hh"
#include "WCSimGeometry.hh"
#include "WCSimInterface.hh"
#include "WCSimLikelihoodTrackBase.hh"
#include "WCSimLikelihoodTrackFactory.hh"
#include "WCSimPiZeroSeed.hh"

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TPrincipal.h>
#include <TVector3.h>
#include <TVectorD.h>

#include <vector>

#ifndef REFLEX_DICTIONARY
ClassImp(WCSimPiZeroElectronAdjuster)
#endif

WCSimPiZeroElectronAdjuster::WCSimPiZeroElectronAdjuster(WCSimFitterConfig * config, WCSimLikelihoodTrackBase * singleElectron, const double &minus2LnL)
 : WCSimPiZeroSeeder(config)
{
	// TODO Auto-generated constructor stub
	fEvent = -999;
	fMadeSeeds = false;
	fSingleElectronTrack = singleElectron;
	fSingleElectron2LnL = minus2LnL;
}

WCSimPiZeroElectronAdjuster::~WCSimPiZeroElectronAdjuster()
{
	// TODO Auto-generated destructor stub
}


void WCSimPiZeroElectronAdjuster::MakeSeeds()
{
	fPiZeroSeeds = IterateOverConversionDistances();
	fMadeSeeds = true;
}


std::vector<WCSimPiZeroSeed*> WCSimPiZeroElectronAdjuster::IterateOverConversionDistances()
{
	std::vector<WCSimPiZeroSeed*> allTracks;

	static const int numDists = 1;
	static const double firstDist = 50.0; //cm
	// static const double secondDist = 250.0; //cm
	double conversionDistancesToTry[numDists] = { firstDist }; //, secondDist };
	for(int iTrack1Dist = 0; iTrack1Dist < numDists; ++iTrack1Dist )
	{
		for(int iTrack2Dist = 0; iTrack2Dist < numDists; ++iTrack2Dist)
		{
			std::vector<WCSimPiZeroSeed* > trackOptions;
			trackOptions = IterateOverFirstTrackPerturbations(conversionDistancesToTry[iTrack1Dist],
															  conversionDistancesToTry[iTrack2Dist]);
			std::cout << "conversionDistances = " << conversionDistancesToTry[iTrack1Dist] << "  " << conversionDistancesToTry[iTrack2Dist] << std::endl;
			std::cout << "size of trackOptions = " << trackOptions.size() << std::endl;
			for(unsigned int iTrack = 0; iTrack < trackOptions.size(); ++iTrack)
			{
				std::cout << "Final seed tracks here are: " << std::endl;
				trackOptions.at(iTrack)->Print();
				allTracks.push_back(trackOptions.at(iTrack));
				std::cout << "Size of allTracks = " << allTracks.size() << std::endl;
			}
		}
	}
	std::cout << "Final size of allTracks = " << allTracks.size() << std::endl;
	return allTracks;
}

std::vector<WCSimPiZeroSeed*> WCSimPiZeroElectronAdjuster::IterateOverFirstTrackPerturbations(const double &convDistTrack1, const double &convDistTrack2)
{
	std::vector<WCSimPiZeroSeed*> allTracks;

	std::vector<TVector3> directionsToTry = GetFirstTrackDirectionsToTry();
	std::vector<TVector3>::iterator directionItr = directionsToTry.begin();

	double currentBestSimilar2LnL = -999.9;
	double currentBestDifferent2LnL = -999.9;
	WCSimPiZeroSeed * bestSimilarTracks = NULL;
	WCSimPiZeroSeed * bestDifferentTracks = NULL;

	while(directionItr != directionsToTry.end())
	{
		std::vector<WCSimPiZeroSeed *> tracks;

		std::cout << "I think the distance is " << std::distance(directionsToTry.begin(), directionItr) <<std::endl;

		tracks = GridSearchOverSecondTrackDirection(
				convDistTrack1, convDistTrack2,
				*directionItr,
				std::distance(directionsToTry.begin(), directionItr));
		std::cout << "Size of vector returned by grid search = " << tracks.size() << std::endl;

		for(unsigned int iSeed = 0; iSeed < tracks.size(); ++iSeed)
		{
			std::cout << "Similar energy? " << SimilarEnergy(tracks.at(iSeed)) << " and minus2LnL = " << tracks.at(iSeed)->GetMinus2LnL() << std::endl;

			if(SimilarEnergy(tracks.at(iSeed)))
			{
				if( tracks.at(iSeed)->GetMinus2LnL() < currentBestSimilar2LnL || (currentBestSimilar2LnL == -999.9 && tracks.at(iSeed)->GetMinus2LnL() != -999.9) )
				{
					std::cout << "Setting bestSimilarTracks" << std::endl;
					bestSimilarTracks = tracks.at(iSeed);
					currentBestSimilar2LnL = bestSimilarTracks->GetMinus2LnL();
				}
				else
				{
					delete tracks.at(iSeed);
					tracks.at(iSeed) = 0x0;
				}
			}
			else
			{
				if( tracks.at(iSeed)->GetMinus2LnL() < currentBestDifferent2LnL || (currentBestDifferent2LnL == -999.9  && tracks.at(iSeed)->GetMinus2LnL() != -999.9) )
				{
					std::cout << "Setting bestDifferentTracks" << std::endl;
					bestDifferentTracks = tracks.at(iSeed);
					currentBestDifferent2LnL = bestDifferentTracks->GetMinus2LnL();
				}
				else
				{
					delete tracks.at(iSeed);
					tracks.at(iSeed) = 0x0;
				}
			}
		}
		tracks.clear();
		++directionItr;
	}

	std::cout << "Is bestSimilarTracks null?" << (bestSimilarTracks == NULL) << std::endl;
	if(bestSimilarTracks != NULL)
	{
		std::cout << "bestSimilarTracks has -2LnL = " << bestSimilarTracks->GetMinus2LnL() << std::endl;
		allTracks.push_back(bestSimilarTracks);
	}
	std::cout << "Is bestDifferentTracks null?" << (bestDifferentTracks == NULL) << std::endl;
	if(bestDifferentTracks != NULL)
	{
		std::cout << "bestDifferentTracks has -2LnL = " << bestDifferentTracks->GetMinus2LnL() << std::endl;
		allTracks.push_back(bestDifferentTracks);
	}


	return allTracks;
}

std::vector<TVector3> WCSimPiZeroElectronAdjuster::GetFirstTrackDirectionsToTry()
{
	// Use the single track seed to define the plane containing its vertex, normal to its direction:
	TVector3 planeContains = fSingleElectronTrack->GetVtx();
	TVector3 planeUnitNormal = fSingleElectronTrack->GetDir().Unit();
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



	// We can't just take the component of the PMT position relative to the track vertex that
	// is perpendicular to the track direction because our ring might straddle an edge, e.g. between
	// the barrel and endcap.  For example, with a track directed along x, everything on the roof
	// would have the same projection
	// Instead we need to project all the PMTs into a plane the same distance away from the vertex
	// in the direction parallel to the track momentum.  Because we are ultimately only interested in the
	// direction, this can be any distance so long as it's constant - I've picked the cylinder radius
	// because it's a convenient scale to visualise
	double projectToDistance = WCSimGeometry::Instance()->GetCylRadius();
	std::vector<TVector2> projections;

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
		myPrincipal.AddRow(point);
		meanX += point[0]/projections.size();
		meanY += point[1]/projections.size();
	}
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
	TVector3 majorAxis = eigenvector.X() * planeAxis1 + eigenvector.Y() * planeAxis2;
	eigenvector = TVector2(eigenvectors(!secondVectorIsMajorAxis,0), eigenvectors(!secondVectorIsMajorAxis,1));
	TVector3 minorAxis = eigenvector.X() * planeAxis1 + eigenvector.Y() * planeAxis2;

	// Normalise them
	majorAxis = majorAxis.Unit();
	minorAxis = minorAxis.Unit();

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
		TVector3 tmp = fSingleElectronTrack->GetDir(); // TVector3::Rotate overwrites the vector, so need a copy
		tmp.Rotate(rotationsAboutMajor[i], majorAxis);
		directionsToTry.push_back(tmp);
	}

	static const int numRotMin = 2; // Rotations about the minor axis
	double rotationsAboutMinor[numRotMin] = { 0.05*M_PI, -0.05 * M_PI };
	for(int jRot = 0; jRot < numRotMin; ++jRot)
	{
		TVector3 tmp = fSingleElectronTrack->GetDir();
		tmp.Rotate(rotationsAboutMinor[jRot], minorAxis);
		directionsToTry.push_back(tmp);
	}

	directionsToTry.clear();
	TVector3 tmp = fSingleElectronTrack->GetDir(); // TVector3::Rotate overwrites the vector, so need a copy
	tmp.Rotate(rotationsAboutMajor[0], majorAxis);
	directionsToTry.push_back(tmp);

	return directionsToTry;
}

std::vector<WCSimPiZeroSeed*> WCSimPiZeroElectronAdjuster::GridSearchOverSecondTrackDirection(
		const double &convDistTrack1,
		const double &convDistTrack2,
		const TVector3 &track1Dir,
		int numTimesAlready)
{
	// std::cout << "*** WCSimPiZeroFitter::GridSearchOverSecondTrackDirection *** " << std::endl;

	// Best log likelihoods for tracks with similar (Ebig / Esmall < 5) and different energies
	double bestSimilar2LnL = -999.9;
	double bestDifferent2LnL = -999.9;

	const unsigned int nBinsTh = 5; // Number of bins in theta
	const unsigned int nBinsPhi = 5; // Number of bins in phi

	WCSimFitterConfig * backupConfig = fFitterConfig;
	WCSimFitterConfig config;
	config.SetNumTracks(2);

	config.SetTrackType(0, "PhotonLike");
	config.SetTrackType(1, "PhotonLike");
	config.SetJoinParametersTogether(0, 1, "kVtxX");
	config.SetJoinParametersTogether(0, 1, "kVtxY");
	config.SetJoinParametersTogether(0, 1, "kVtxZ");
	config.SetJoinParametersTogether(0, 1, "kVtxT");
	if(fFitterConfig->GetIsPiZeroFit())
	{
		config.SetIsPiZeroFit(fFitterConfig->GetIsPiZeroFit());
		config.SetForcePiZeroMass(fFitterConfig->GetForcePiZeroMass());
	}
	config.SetParameter(0, "kVtxX", fFitterConfig->GetParMin(0, "kVtxX"), fFitterConfig->GetParMax(0, "kVtxX"), fSingleElectronTrack->GetX(), false);
	config.SetParameter(0, "kVtxY", fFitterConfig->GetParMin(0, "kVtxY"), fFitterConfig->GetParMax(0, "kVtxY"), fSingleElectronTrack->GetY(), false);
	config.SetParameter(0, "kVtxZ", fFitterConfig->GetParMin(0, "kVtxZ"), fFitterConfig->GetParMax(0, "kVtxZ"), fSingleElectronTrack->GetZ(), false);
	config.SetParameter(0, "kVtxT", fFitterConfig->GetParMin(0, "kVtxT"), fFitterConfig->GetParMax(0, "kVtxT"), fSingleElectronTrack->GetT(), false);
	config.SetParameter(0, "kDirTh", fFitterConfig->GetParMin(0, "kDirTh"), fFitterConfig->GetParMax(0, "kDirTh"), track1Dir.Theta(), false);
	config.SetParameter(0, "kDirPhi", fFitterConfig->GetParMin(0, "kDirPhi"), fFitterConfig->GetParMax(0, "kDirPhi"), track1Dir.Phi(), false);
	config.SetParameter(0, "kEnergy", fFitterConfig->GetParMin(0, "kEnergy"), fFitterConfig->GetParMax(0, "kEnergy"), fSingleElectronTrack->GetE() , false);
	config.SetParameter(0, "kConversionDistance", fFitterConfig->GetParMin(0, "kConversionDistance"), fFitterConfig->GetParMax(0, "kConversionDistance"), convDistTrack1, false);
	config.SetParameter(1, "kDirTh", fFitterConfig->GetParMin(1, "kDirTh"), fFitterConfig->GetParMax(1, "kDirTh"), fFitterConfig->GetParStart(1, "kDirTh"), false);
	config.SetParameter(1, "kDirPhi", fFitterConfig->GetParMin(1, "kDirPhi"), fFitterConfig->GetParMax(1, "kDirPhi"), fFitterConfig->GetParStart(1, "kDirPhi"), false);
	config.SetParameter(1, "kEnergy", fFitterConfig->GetParMin(1, "kEnergy"), fFitterConfig->GetParMax(1, "kEnergy"), fFitterConfig->GetParStart(1, "kEnergy"), false);
	config.SetParameter(1, "kConversionDistance", fFitterConfig->GetParMin(1, "kConversionDistance"), fFitterConfig->GetParMax(1, "kConversionDistance"), convDistTrack2, false);

	fFitterConfig = &config;
	fFitterTrackParMap = WCSimFitterTrackParMap(&config);
	fFitterTrackParMap.Set();

	// Default the best values to something junk so we can tell when they've been set
	double bestSimilarTheta = -999.9;
	double bestSimilarPhi = -999.9;
	double bestSimilarEnergy1 = -999.9;
	double bestSimilarEnergy2 = -999.9;

	double bestDifferentTheta = -999.9;
	double bestDifferentPhi = -999.9;
	double bestDifferentEnergy1 = -999.9;
	double bestDifferentEnergy2 = -999.9;

	// Vector orthogonal to the direction of track 1: we do the theta rotation about this axis
	// and the phi rotation about track1Dir
	TVector3 orthog = track1Dir.Orthogonal();

	// Start the grid search
	for(unsigned int thetaBin = 1; thetaBin <= nBinsTh; ++thetaBin)
	{
		double dTheta = (0.125 * M_PI * thetaBin)/nBinsTh;


		for(unsigned int phiBin = 1; phiBin <= nBinsPhi; ++phiBin)
		{
			double dPhi = phiBin * 2.0 * M_PI/nBinsPhi;
			TVector3 track2Dir = track1Dir;
			track2Dir.Rotate(dTheta, orthog);
			track2Dir.Rotate(dPhi, track1Dir);

			std::cout << "Grid search number " << numTimesAlready+1 << "  Bin " << phiBin + (thetaBin-1) * nBinsPhi << "/" << nBinsTh*nBinsPhi << std::endl;
			std::cout << "Angle between two tracks = " << track2Dir.Angle(track1Dir) << "   (" << track2Dir.Angle(track1Dir) * 180.0 / M_PI << ") deg" << std::endl;

			fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirTh, track2Dir.Theta());
			fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kDirPhi, track2Dir.Phi());

			std::pair<double, double> energies = GetPiZeroPhotonEnergies(track1Dir, fSingleElectronTrack->GetE(), track2Dir);
			if(energies.second > 5000)
			{
				continue;
			}

			if(   numTimesAlready > 0
			   && ( energies.second < 0.05 * energies.first || energies.first < 0.05 * energies.second ))
			{
				continue;
			}
			fFitterTrackParMap.SetCurrentValue(0, FitterParameterType::kEnergy, energies.first);
			fFitterTrackParMap.SetCurrentValue(1, FitterParameterType::kEnergy, energies.second);

			std::vector<double> startVals = fFitterTrackParMap.GetCurrentValues();
			Double_t * x = &(startVals[0]);
     		double tempMin = WrapFunc(x);

     		if(!SimilarEnergy(energies.first, energies.second))
     		{
     			if( tempMin < bestDifferent2LnL || (bestDifferent2LnL < 0 && tempMin > 0) )
     			{
     				std::cout << "Found new best different lnl of " << tempMin << std::endl;
     				bestDifferent2LnL = tempMin;
     				bestDifferentTheta = track2Dir.Theta();
     				bestDifferentPhi = track2Dir.Phi();
     				bestDifferentEnergy1 = energies.first;
     				bestDifferentEnergy2 = energies.second;
     			}
     		}
     		else
     		{
     			if( tempMin < bestSimilar2LnL || (bestSimilar2LnL < 0 && tempMin > 0) )
				{
     				std::cout << "Found new best similar lnl of " << tempMin << std::endl;
     				bestSimilar2LnL = tempMin;
     				bestSimilarTheta = track2Dir.Theta();
     				bestSimilarPhi = track2Dir.Phi();
     				bestSimilarEnergy1 = energies.first;
     				bestSimilarEnergy2 = energies.second;
				}
     		}
		}
	}

	// Return the best-fit
	std::vector<WCSimPiZeroSeed*> best;
	if( bestSimilar2LnL > 0)
	{
		std::map<FitterParameterType::Type, double> extraPars;
		extraPars[FitterParameterType::kConversionDistance] = convDistTrack1;
		WCSimLikelihoodTrackBase * firstTrack = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::PhotonLike,
				fSingleElectronTrack->GetX(), fSingleElectronTrack->GetY(), fSingleElectronTrack->GetZ(),
				fSingleElectronTrack->GetT(),
				track1Dir.Theta(), track1Dir.Phi(),
				bestSimilarEnergy1,
				extraPars);


		std::map<FitterParameterType::Type, double> extraPars2;
		extraPars2[FitterParameterType::kConversionDistance] = convDistTrack2;
		WCSimLikelihoodTrackBase * secondTrack = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::PhotonLike,
				fSingleElectronTrack->GetX(), fSingleElectronTrack->GetY(), fSingleElectronTrack->GetZ(),
				fSingleElectronTrack->GetT(),
				bestSimilarTheta,
				bestSimilarPhi,
				bestSimilarEnergy2,
				extraPars2);
		WCSimPiZeroSeed * seed = new WCSimPiZeroSeed(firstTrack, secondTrack, bestSimilar2LnL);
		best.push_back(seed);
		delete firstTrack;
		delete secondTrack;
	}

  if( bestDifferent2LnL > 0)
  {
		std::map<FitterParameterType::Type, double> extraPars;
		extraPars[FitterParameterType::kConversionDistance] = convDistTrack1;
		WCSimLikelihoodTrackBase * firstTrack = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::PhotonLike,
				fSingleElectronTrack->GetX(), fSingleElectronTrack->GetY(), fSingleElectronTrack->GetZ(),
				fSingleElectronTrack->GetT(),
				track1Dir.Theta(), track1Dir.Phi(),
				bestDifferentEnergy1,
				extraPars);


		std::map<FitterParameterType::Type, double> extraPars2;
		extraPars2[FitterParameterType::kConversionDistance] = convDistTrack2;
		WCSimLikelihoodTrackBase * secondTrack = WCSimLikelihoodTrackFactory::MakeTrack(
				TrackType::PhotonLike,
				fSingleElectronTrack->GetX(), fSingleElectronTrack->GetY(), fSingleElectronTrack->GetZ(),
				fSingleElectronTrack->GetT(),
				bestDifferentTheta,
				bestDifferentPhi,
				bestDifferentEnergy2,
				extraPars2);
		WCSimPiZeroSeed * seed = new WCSimPiZeroSeed(firstTrack, secondTrack, bestSimilar2LnL);
		best.push_back(seed);
		delete firstTrack;
		delete secondTrack;
  }
	std::cout << "Best similar 2LnL was " << bestSimilar2LnL << std::endl;
	std::cout << "Best different 2LnL was " << bestDifferent2LnL << std::endl;
	std::cout << "Size of best-fit vector was " << best.size() << std::endl;
	std::cout << "Printing bestfit tracks:" << std::endl;

	for(unsigned int i = 0; i < best.size(); ++i)
	{
		std::cout << "Track pair " << i << std::endl;
		best.at(i)->Print();
	}

	fFitterConfig = backupConfig;
	return best;
}


// Get the energy of photon track 1 given the energy of photon track 0
// assuming they arise from the decay of a pi zero (and so have its invariant
// mass)
// \param x The array of current values that will be passed to/returned from WrapFuncPiZero
//          (which comes from the FitterTrackParMap)
std::pair<double, double> WCSimPiZeroElectronAdjuster::GetPiZeroPhotonEnergies(const TVector3 &track1Dir, const double &track1Energy, const TVector3 &track2Dir )
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

  double cosTheta12 =   sinThTrk1 * cosPhiTrk1 * sinThTrk2 * cosPhiTrk2
                      + sinThTrk1 * sinPhiTrk1 * sinThTrk2 * sinPhiTrk2
                      + cosThTrk1 * cosThTrk2;
                      // This is the dot product of unit vectors in the directions of track 0 and track 1
                      // i.e. the cosine of the angle between the two photons

  double track1Estimator = GetFirstTrackEnergyEstimator(track1Energy, cosTheta12);
  double track2Estimator = 0.5 * massOfPiZeroMeV * massOfPiZeroMeV / (track1Estimator * (1.0 - cosTheta12));
  /*if(track1Estimator < 100 && track2Estimator > 100)
  {
    double scale = 100.0 / track1Estimator;
    track1Estimator /= scale;
    track2Estimator *= scale;
  }
  else if(track2Estimator < 100 && track1Estimator > 100)
  {
    double scale = 100.0 / track1Estimator;
    track1Estimator *= scale;
    track2Estimator /= scale;
  }*/

  TLorentzVector fourMom0( track1Estimator * sinThTrk1 * cosPhiTrk1, track1Estimator * sinThTrk1 * sinPhiTrk1, track1Estimator * cosThTrk1, track1Estimator);
  TLorentzVector fourMom1( track2Estimator * sinThTrk2 * cosPhiTrk2, track2Estimator * sinThTrk2 * sinPhiTrk2, track2Estimator * cosThTrk2, track2Estimator);
  // std::cout << "E1 = " << track1Energy << "   E2 = " << track2Energy << "   Separation = " << TMath::ACos(cosTheta12) << "  Total four-momentum = " << (fourMom0 + fourMom1).Mag() << std::endl;
  std::cout << "Fit energy = " << track1Energy << "   so set E1 = " << track1Estimator << " and E2 = " << track2Estimator << std::endl;
  std::cout << "Separation = " << acos(cosTheta12) << "  so four-momentum = " << (fourMom0 + fourMom1).Mag() << std::endl;

  return std::make_pair(track1Estimator, track2Estimator);
}

double WCSimPiZeroElectronAdjuster::GetFirstTrackEnergyEstimator(const double &singleTrackFitE, const double &cosThetaSep)
{
  double theta = acos(cosThetaSep);
  while(theta < 0){
    theta += 2*M_PI;
  }
  double rmtor = 0.549803 - 0.148013 * theta; // (Reco - true) over reco, from an empirical fit
  double energy = singleTrackFitE * (1 - rmtor);
  if(energy < 50)
  {
    energy = 50;
  }
  return energy;
}

bool WCSimPiZeroElectronAdjuster::SimilarEnergy(WCSimPiZeroSeed* seed) const
{
	return SimilarEnergy(seed->GetTrack1(), seed->GetTrack2());
}

bool WCSimPiZeroElectronAdjuster::SimilarEnergy(
		WCSimLikelihoodTrackBase* track1, WCSimLikelihoodTrackBase* track2) const
{
	return SimilarEnergy(track1->GetE(), track2->GetE());
}

bool WCSimPiZeroElectronAdjuster::SimilarEnergy(const double &energy1, const double &energy2) const
{
	float factor = 2.0;
	return ( energy1 < (factor * energy2) || energy2 < (factor * energy1) );
}
