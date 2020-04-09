#include "TCanvas.h"
#include "TDatime.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
double IntegrateMeanTimeToMeanFirstTime(double * x, double * par);
double IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons);
double ConvertMeanTimeToMeanFirstTime(double * x, double * par);
double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons);
TH1D * GetMCHist(double mean, double sigma, int nPhotons, int nTrials);

void wc_plotFirstTimePDF(){
  double n      = 3;
  double sigma  = 1.0;
  double sqrtPi = TMath::Sqrt(TMath::Pi());
  double mean   = 100;
  int nTimes    = 10000;

  // MC histogram of the first arrival time when n photons arrive with a given mean and rms arrival time
  TH1D * hMonteCarlo = GetMCHist(mean, sigma, n, nTimes);
  hMonteCarlo->GetXaxis()->CenterTitle();
  hMonteCarlo->GetYaxis()->CenterTitle();
  std::cout << hMonteCarlo->Integral("width") << std::endl;


  // The PDF for the first arrival time
  TF1 * fArrival = new TF1("fArrival", ConvertMeanTimeToMeanFirstTime, mean - 5 * sigma, mean + 5 * sigma, 3);
  fArrival->SetParameters(n, mean, sigma);

  // The integral of the PDF for the first arrival time
  TF1 * fInt = new TF1("fInt", IntegrateMeanTimeToMeanFirstTime, mean - 5 * sigma, mean + 5 * sigma, 3);
  fInt->SetParameters(n, mean, sigma);

  // Format and then draw the histogram and the two functions
  TCanvas * can = new TCanvas("can","",800, 600);
  hMonteCarlo->SetLineColor(kAzure);
  hMonteCarlo->SetFillColor(kAzure);
  hMonteCarlo->SetFillStyle(3002);
  hMonteCarlo->SetStats(0);
  hMonteCarlo->SetMaximum(1.2);
  fArrival->SetLineColor(kRed);
  fInt->SetLineColor(kViolet-1);
  hMonteCarlo->Draw();
  fInt->Draw("SAME");
  fArrival->Draw("SAME");

  // Build legend
  TLegend * leg = new TLegend(0.22,0.85,0.47,0.6);
  leg->AddEntry(hMonteCarlo, "Toy MC", "LF");
  leg->AddEntry(fArrival, "PDF of first arrival time", "L");
  leg->AddEntry(fInt, "CDF of first arrival time", "L");
  leg->SetTextFont(42);
  leg->Draw();


  // Save in useful formats
  can->SaveAs("predictedFirstArrivalTime.png");
  can->SaveAs("predictedFirstArrivalTime.pdf");
  can->SaveAs("predictedFirstArrivalTime.root");
  can->SaveAs("predictedFirstArrivalTime.C");

  return;
}



/**
 * @brief Build toy MC histogram for the distribution of the first arrival time of photons
 *
 * @param mean Mean arrival time of all the photons
 * @param sigma Gaussian width of photon arrival times
 * @param nPhotons Total number of photons
 * @param nTrials Number of times to run the MC
 *
 * @return Histogram of the earliest of the n arrival times, scaled to unit area
 */
TH1D * GetMCHist(double mean, double sigma, int nPhotons, int nTrials)
{
	TRandom3 rndm;
	TH1D * hist = new TH1D("hMonteCarlo",";Time (ns);Number of events", 100, mean-5*sigma, mean+5*sigma);

	for(int trial = 0; trial < nTrials; ++trial)
	{
		double min = -999.9;
		for(int i = 0 ; i < nPhotons; ++i)
		{
			double ran = rndm.Gaus(mean, sigma);
			if(ran < min || i == 0)
			{
				min = ran;
			}
		}
		hist->Fill(min);
	}
	hist->Scale(1.0/hist->Integral("width"));

	return hist;

}

/**
 * Given n samples drawn from a Gaussian distribution, return the PDF for the earliest sample,
 * with double pointer arguments for making it into a TF1
 *
 * @param x Pointer to the value at which we want to evaluate the PDF
 * @param par Array containing the underlying Gaussian's mean and RMS, and number of photons
 * @return The PDF for the first arrival time being equal to x[0]
 * @param par Array containing the underlying Gaussian's mean and RMS, and number of photons
 */
double ConvertMeanTimeToMeanFirstTime(double * x, double * par)
{
	double nPhotons = par[0];
	double mean = par[1];
	double sigma = par[2];
	double t = x[0];
	return ConvertMeanTimeToMeanFirstTime(t, mean, sigma, nPhotons);
}

/**
 * Evaluate the integral of the PDF for the first arrival time
 * @param t Integral evaluated from -infinity to t
 * @param mean Mean of the underlying Gaussian
 * @param sigma Width of the underlying Gaussian
 * @param nPhotons Number of photons, i.e. number of times that Gaussian is sampled
 * @return Integral of the PDF for the first arrival time from -infinity to t
 */
double IntegrateMeanTimeToMeanFirstTime( const double &t, const double &mean, const double &sigma, const double &nPhotons)
{
  double integral(0.0);
  if(nPhotons <= 25)
  {
    integral = 1.0 - 1.0/(TMath::Power(2,nPhotons)) * TMath::Power((1.0 - TMath::Erf((t - mean)/(TMath::Sqrt2() * sigma))), nPhotons);
  }
  else
  {
    double exponent = nPhotons * (TMath::Log(1 - TMath::Erf((t - mean) / (TMath::Sqrt2() * sigma))) - TMath::Log(2));
    integral = 1.0 - TMath::Exp(exponent);
  }
  return integral;
}

/**
 * Evaluate the integral of the PDF for the first arrival time - with double
 * pointer arguments so it can be made into a TF1
 * @param x Integral evaluated from -infinity to x[0]
 * @param par Array containing the underlying Gaussian's mean and RMS, and number of photons
 * @return Integral of the PDF for the first arrival time from -infinity to x[0]
 */
double IntegrateMeanTimeToMeanFirstTime(double * x, double * par)
{
  double &nPhotons = (par[0]);
  double &mean = (par[1]);
  double &sigma = (par[2]);
  double &t = (x[0]);

  return IntegrateMeanTimeToMeanFirstTime(t, mean, sigma, nPhotons);
}

/**
 * Given n samples drawn from a Gaussian distribution, return the PDF for the earliest sample,
 * with double pointer arguments for making it into a TF1
 *
 * @param t Time at which we want to evaluate the PDF
 * @param mean Mean arrival time of the photons
 * @param sigma Gaussian width of photon arrival times
 * @param nPhotons Total number of photons
 * @return The PDF for the first arrival time being t
 */
double ConvertMeanTimeToMeanFirstTime(const double &t, const double &mean, const double &sigma, const double &nPhotons)
{
	std::cout << "t = " << t << " mean  = " << mean << "  sigma = " << sigma << " nPhotons = " << nPhotons << std::endl;
    double prob = 0.0;
    if(nPhotons == 0 || sigma == 0) { return 0; }
    if(nPhotons < 25)
    {
	    prob =   nPhotons * TMath::Sqrt2() / (TMath::Power(2, nPhotons) * sigma * TMath::Sqrt(TMath::Pi()))
		 	  * TMath::Power(1.0 - TMath::Erf((t - mean)/(TMath::Sqrt2()*sigma)), nPhotons-1)
		   	  * TMath::Exp(-(t - mean)*(t - mean) / (2 * sigma * sigma));
    }
    else // Use logarithms to take care of the two large n powers almost cancelling
    {
        double logL =  TMath::Log(nPhotons * TMath::Sqrt2() / (sigma * TMath::Sqrt(TMath::Pi()) ))
                     - nPhotons * TMath::Log(2.)
                     + (nPhotons - 1) * TMath::Log(1.0 - TMath::Erf((t - mean) / (2 * sigma * sigma)))
                     - (t - mean) * (t - mean) / (2 * sigma * sigma);
        prob = TMath::Exp(logL);

    }
	return prob;
}
