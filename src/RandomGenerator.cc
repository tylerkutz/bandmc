#include "RandomGenerator.hh"

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"

#include "constants.hh"
#include "GenTree.hh"

using namespace std;


RandomGenerator::RandomGenerator() {

	minTrueKE = 0.015;
	maxTrueKE = 1.; 
	minMomR = 0.2;
	maxMomR = 0.65;

	maxThetaR = 178.*M_PI/180.;
	minThetaR = 150.*M_PI/180.;


	bandZ = -267; // cm

	fRand = new TRandom3(0);
	
	maxBetaR = beta(maxMomR);
	minBetaR = beta(minMomR);
	minCosThetaR = cos(maxThetaR);
	maxCosThetaR = cos(minThetaR);
	minTime = -100;//-bandZ/(maxBetaR*cAir) - 5.;  		// Put in 5 ns for cushion
	maxTime = 100;//(7.2*5-bandZ)/(minBetaR*cAir) + 5.; 	// Put in another 5 ns for more cushion
	timeWindow = maxTime - minTime;
	BAND_center =  bandZ - (7.2*5)/2.;

	inclFile = new TFile("dat/inclusive_sample_background.root");
	inclTree = (TTree*)inclFile->Get("T");

	inclTree->SetBranchAddress("p_e", 	&p_e);
	inclTree->SetBranchAddress("theta_e", 	&theta_e);
	inclTree->SetBranchAddress("phi_e",	&phi_e);

	nBGSamp = inclTree->GetEntries();

	nBG = 0; 

	step = fRand->Integer(100);

}

void RandomGenerator::GenerateRandom(Gen_Event* thisEvent) {

	// Generate random neutron
	double cosThetaR = minCosThetaR + fRand->Rndm()*(maxCosThetaR-minCosThetaR);
	double phiR = 2.*M_PI * fRand->Rndm();
	double keR = getNeutronKE(fRand->Rndm());

	// Calculate the derived neutron info
	double thetaR = acos(cosThetaR);
	double momR = sqrt(sq(keR + mN) - mN*mN);
	double betaR = momR / (keR + mN);

	thisEvent->particles.clear();
	
	// Need to get electron info from inclusive skim
	inclTree->GetEntry(nBG);
	nBG += step;
	nBG = nBG%nBGSamp;

	// Electron info
	Gen_Particle electron;
	electron.type="electron";
	electron.momentum.SetMagThetaPhi(p_e, theta_e, phi_e);
	thisEvent->particles.push_back(electron);

	// Neutron info
	Gen_Particle neutron;
	neutron.type="neutron";
	neutron.momentum.SetMagThetaPhi(momR,thetaR,phiR);
	thisEvent->particles.push_back(neutron);

	/*
	// Update cross section info
	double neutronCS = getNeutronCS(minTrueKE) * (2.*M_PI)*(maxCosThetaR - minCosThetaR);
	TVectorT<double> csSqVec(3);
	csSqVec[0]=neutronCS*timeWindow*(*csVec)[0]; // Units of nb^2 * ns
	csSqVec[1]=neutronCS*timeWindow*(*csVec)[1];
	csSqVec[2]=timeWindow;
	csSqVec.Write("totalCSSq");
	*/

}

double RandomGenerator::getNeutronCS(double threshold) {

	double pavel_rate = 967451. * exp(-29.2539 * threshold) + 5.43632e+06 * exp(-548.598 * threshold);

	return pavel_rate / (3.E36) * (1.E33) / (0.1); // divide by Pavel's lumi, convert to nb, divide by 0.1 sr
}

double RandomGenerator::CDF(double Tr)
{
	return 1. - getNeutronCS(Tr)/getNeutronCS(minTrueKE);
}

double RandomGenerator::getNeutronKE(double r)
{
	double minKE = minTrueKE;
	double maxKE = maxTrueKE;
	double testKE = 0.5*(minKE + maxKE);

	while (fabs(maxKE - minKE) > 1.E-6)
	{
		double testCDF = CDF(testKE);

		if (r < testCDF)
			maxKE = testKE;
		else
			minKE = testKE;

		testKE = 0.5*(minKE + maxKE);
	}

	return testKE;
}

double RandomGenerator::sq(double x){ 
	return x*x; 
}

double RandomGenerator::beta(double p){ 
	return 1./sqrt(1. + sq(mN/p)); 
};

