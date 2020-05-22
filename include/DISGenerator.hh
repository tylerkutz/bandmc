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

/*
	double sigmainput;          //sigma parameter in rescattering amplitude [mb], default value
	double betainput;            //beta parameter in rescattering amplitude [GeV^-2], default value
	double epsinput;           //epsilon parameter in rescattering amplitude [], default value
	double lambdainput;
	double betaoffinput;
	int offshellset;              // 0=mass diff suppression, 1=dipole suppression, 2=beta parameterization, 3=no offshell, 4=full off-shell
	int symm;                     //symmetric FSI in inclusive reaction
	int phiavg;                   //average cross section over phi
	int F_param;                  //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

	// Foam ranges
	double min_theta_e; 
	double max_theta_e;
	double min_p_e; 	   
	double max_p_e; 	   

	double min_theta_r; 
	double max_theta_r;
	double min_p_r; 	   
	double max_p_r; 	   

	double min_phi_r;   
	double max_phi_r;   
*/
};


class DISGenerator {

public:
	DISGenerator();
	~DISGenerator();

	void GenerateEvent(Gen_Event*);
	void GetTotalCS(double*);

private:

//	void SetDIS();


/*
	// Settings for cross section: 
	double sigmainput;          //sigma parameter in rescattering amplitude [mb], default value
	double betainput;            //beta parameter in rescattering amplitude [GeV^-2], default value
	double epsinput;           //epsilon parameter in rescattering amplitude [], default value
	double lambdainput;
	double betaoffinput;
	int offshellset;              // 0=mass diff suppression, 1=dipole suppression, 2=beta parameterization, 3=no offshell, 4=full off-shell
	int symm;                     //symmetric FSI in inclusive reaction
	int phiavg;                   //average cross section over phi
	int F_param;                  //structure functions parametrization 0=SLAC, 1=Christy&Bosted, 2=Alekhin et al leading twist

	// Foam ranges
	double min_theta_e; 
	double max_theta_e;
	double min_p_e; 	   
	double max_p_e; 	   

	double min_theta_r; 
	double max_theta_r;
	double min_p_r; 	   
	double max_p_r; 	   

	double min_phi_r;   
	double max_phi_r;   
*/
	disCS* fDIS;
	TRandom3* rand;
	TFoam* csFoam;

};

#endif


