#ifndef SIMULATE_HH
#define SIMULATE_HH

#include "TLorentzVector.h"

class TVector3;
class TRandom3;

class RadGen;
class Smearing;
class IO;
class ElectronFiducial;
class BAND;

class Simulate{

public:
	Simulate(int, int);
	~Simulate();
	int SimulateEvent(TVector3, TVector3); 
	int SimulateBackground(TVector3, TVector3); 
	void SetIO(IO* io) {fIO = io;}

private:

	int nSim;
	int nBG;

	int doRadiation;
	int doSmearing;

	double targetL;
	double barThickness;
	double zUp;
	double zDown;

	RadGen* fRad;
	Smearing* fSmear;
	IO* fIO;	
	ElectronFiducial* fFiducial;
	BAND* fBAND;

	TRandom3* fRand;

	void SetData();

	double Me;
	double Mp;
	double Mn;
	double MD;
	double eBeam;
	double c;

	double E_e;

	TLorentzVector e4;

	double p_e; 		
	double theta_e; 	
	double phi_e;		
	double p_n;		
	double theta_n;	
	double phi_n;		
	double E_n;
	double Q2;		
	double Q2rad;		
	double Q2true;		
	double xB;		
	double W2;		
	double Ebeam;		
	double gated_charge;	
	double livetime;	
	double starttime;	
	double ePid;		
	double eCharge;	
	double eStatus;	
	double eTime;		
	double eBeta;		
	double eChi2pid;	
	double E_tot;		
	double E_pcal;	
	double t_e;		
	double dL_e;		
	double lU;		
	double lV;		
	double lW;		
	double e_vtx;		
	double e_vty;		
	double e_vtz;		
	double q;		
	double theta_q;	
	double phi_q;		
	double nu;		
	double nMult;		
	double barID;		
	double dL_n;		
	double nTime;		
	double nBeta;
	double nEdep;		
	int bg;
	double CosTheta_nq;
	double Xp;
	double Wp;
	double As;

	double radweight;

};

#endif

