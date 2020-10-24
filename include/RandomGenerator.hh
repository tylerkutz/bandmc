#ifndef RANDOMGENERATOR_HH
#define RANDOMGENERATOR_HH

class Gen_Event;
class Gen_Particle;
class TRandom3;
class TFile;
class TTree;

class RandomGenerator{

public:

	RandomGenerator();
	~RandomGenerator();

	void GenerateRandom(Gen_Event*);

private:

	double minTrueKE;
	double maxTrueKE;
	double minMomR;
	double maxMomR;
	double maxThetaR;
	double minThetaR;
	double bandZ; // cm

	double maxBetaR;
	double minBetaR;
	double minCosThetaR;
	double maxCosThetaR;
	double minTime;
	double maxTime;
	double timeWindow;
	double BAND_center;

	double sq(double);
	double beta(double);
	double getNeutronCS(double); // in nb/sr
	double getNeutronKE(double);
	double CDF(double);

	TRandom3* fRand;
	TFile* inclFile;
	TTree* inclTree;

	double pe[3];
//	double p_e;
//	double theta_e;
//	double phi_e;

	int nBGSamp;
	int nBG;
	int step;	


};

#endif

