#ifndef BAND_HH
#define BAND_HH

class TRandom3;
class TF1;
class TH1F;
class TFile;
class TVector3;

class BAND{

public:
	BAND();
	~BAND();
	double GetNeutronAcceptance(double, double, double, double);
	double PointsToBAND(double, double, double);
	double Efficiency(double);
	double GetEdep(double);
	void Smear(TVector3*);

private:
	TRandom3* fRand;	
	TF1* nEffFunc;

	TFile* edepFile;
	double dist_p_n[12];
	double bound_p_n[13];
	TH1F* edep_dist[12];

};

#endif

