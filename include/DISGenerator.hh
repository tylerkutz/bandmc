#ifndef DISGENERATOR_HH
#define DISGENERATOR_HH

#include "TFoamIntegrand.h"

class TRandom3;
class TFoam;
class Gen_Event;
class Radiation;

// Foam Integrand
class disCS : public TFoamIntegrand {

public:
	double Density(int nDim, double * args);
	virtual ~disCS() {};
	inline void SetRadiation(Radiation* rad) {fRad = rad;};
	inline void DoRadiate(int dr) {doRadiation = dr;};

private:
	Radiation* fRad;
	int doRadiation;

};


class DISGenerator {

public:
	DISGenerator(int);
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


