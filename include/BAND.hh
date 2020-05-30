#ifndef BAND_HH
#define BAND_HH

class TRandom3;

class BAND{

public:
	BAND();
	~BAND();
	double GetNeutronAcceptance(double, double, double, double);


private:

	double Efficiency(double);
	double PointsToBAND(double, double, double);

	TRandom3* fRand;	

};

#endif

