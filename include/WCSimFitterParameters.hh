/*
 * WCSimFitterParameters.hh
 *
 *  Created on: 31 Oct 2014
 *      Author: andy
 */

#ifndef WCSIMFITTERPARAMETERS_HH_
#define WCSIMFITTERPARAMETERS_HH_
#include "WCSimLikelihoodTrack.hh"
#include "WCSimTrackParameterEnums.hh"
#include <string>
#include <cassert>
#include <map>
#include <vector>
#include <iostream>

class WCSimFitterParameter {
	public:
		WCSimFitterParameter();

		WCSimFitterParameter(FitterParameterType::Type type, bool isFixed, double start, double min, double max,
				double step = 0.0);

		virtual ~WCSimFitterParameter();

		bool GetIsFixed() const {
			return fIsFixed;
		}

		void SetIsFixed(bool isFixed) {
			fIsFixed = isFixed;
		}

		double GetMax() const {
			return fMax;
		}

		void SetMax(double max) {
			fMax = max;
		}

		double GetMin() const {
			return fMin;
		}

		void SetMin(double min) {
			fMin = min;
		}

		double GetStart() const {
			return fStart;
		}

		void SetStart(double start) {
			fStart = start;
		}

		FitterParameterType::Type GetType() const {
			return fType;
		}

		void SetType(FitterParameterType::Type type) {
			fType = type;
		}

		void SetStep(double step) {
			fStep = step;
		}

		double GetStep() const {
			return fStep;
		}

		void Print();

	private:
		FitterParameterType::Type fType;
		bool fIsFixed;
		double fStart;
		double fMin;
		double fMax;
		double fStep;
};

class WCSimFitterSingleTrackParameters {

	public:
		WCSimFitterSingleTrackParameters();
		virtual ~WCSimFitterSingleTrackParameters();
		void SetDefaultParameters();
		WCSimFitterParameter GetParameter(FitterParameterType::Type type);
		void SetParameter(FitterParameterType::Type type, bool isFixed, double start, double min, double max,
				double step = 0);

		void SetParMin(FitterParameterType::Type type, double min);
		void SetParMax(FitterParameterType::Type type, double max);
		void SetParStart(FitterParameterType::Type type, double start);
		void SetParRange(FitterParameterType::Type type, double min, double max);
		void SetParIsFixed(FitterParameterType::Type type, bool fixIt);
		void SetParStep(FitterParameterType::Type type, double step);

		bool GetParIsFixed(FitterParameterType type);
		double GetParMax(FitterParameterType type);
		double GetParMin(FitterParameterType type);
		double GetParStart(FitterParameterType type);
		double GetParStep(FitterParameterType type);

		unsigned int GetNumParameters();

	private:
		std::map<FitterParameterType::Type, WCSimFitterParameter> fParameters;
};

class WCSimFitterParameters {
	public:
		WCSimFitterParameters();
		virtual ~WCSimFitterParameters();
		void SetNumTracks(unsigned int nTracks);
		void SetTrackType(unsigned int nTrack, TrackType::Type trackType);
		TrackType::Type GetTrackType(const unsigned int &nTrack) const;
		unsigned int GetNumTracks() const {
			return fNumTracks;
		}
		;
		unsigned int GetNumParameters() const;
		unsigned int GetNumIndependentParameters() const;
		void AddTrack(WCSimFitterSingleTrackParameters trackPars);
		WCSimFitterSingleTrackParameters * GetTrackParameters(unsigned int trackNum);

		void JoinParametersTogether(unsigned int track1, unsigned int track2, FitterParameterType::Type type);
		bool GetJoinParametersTogether(unsigned int track1, unsigned int track2, FitterParameterType::Type type);
		bool GetIsParameterJoined(unsigned int track, FitterParameterType::Type type);
		unsigned int GetTrackIsJoinedWith(unsigned int track, FitterParameterType::Type);

	private:
		unsigned int fNumTracks;
		unsigned int fNumParameters;
		std::vector<WCSimFitterSingleTrackParameters> fTrackPars;
		std::map<std::pair<unsigned int, unsigned int>, std::vector<FitterParameterType::Type> > fJoinedParams;
		std::vector<TrackType::Type> fTrackTypes;
};

#endif /* WCSIMFITTERPARAMETERS_HH_ */
