#ifndef DISGENERATOR_HH
#define DISGENERATOR_HH

#include "TFoamIntegrand.h"

class TRandom3;
class TFoam;
class Gen_Event;

// Foam Integrand
class disCS : public TFoamIntegrand {

public:
	double Density(int nDim, double * args);
	virtual ~disCS() {};

};


class DISGenerator {

public:
	DISGenerator();
	~DISGenerator();

	void GenerateEvent(Gen_Event*);
	void GenerateBackground(Gen_Event*);
	void GetTotalCS(double*);

private:

	disCS* fDIS;
	TRandom3* rand;
	TFoam* csFoam;

};

#endif


