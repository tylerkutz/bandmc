#ifndef RADIATION_HH
#define RADIATION_HH

class TGraph2D;
class TString;

class Radiation{

public:
	Radiation();
	~Radiation();

	double GetRadiativeCorrection(double, double);

private:

	void LoadGraph(TString);
	TGraph2D* fRadGraph;

};

#endif

