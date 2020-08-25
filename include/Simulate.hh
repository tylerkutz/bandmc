#ifndef SIMULATE_HH
#define SIMULATE_HH

#include "TLorentzVector.h"

class TVector3;
class TRandom3;

class IO;
class CLAS;
class BAND;

class Simulate{

public:
	Simulate(int, bool);
	~Simulate();
	int SimulateEvent(TVector3, TVector3); 
	int SimulateBackground(TVector3, TVector3); 
	void SetIO(IO* io) {fIO = io;}

private:

	int nSim;
	int nBG;

	bool doSmearing;
	bool isInclusive;

	double targetL;
	double targetR;
	double barThickness;
	double zUp;
	double zDown;

	IO* fIO;	
	CLAS* fCLAS;
	BAND* fBAND;

	TRandom3* fRand;

	void SetEventData();

	double Me;
	double Mp;
	double Mn;
	double MD;
	double eBeam;
	double c;

	double E_e;

	TLorentzVector e4;

	double vx;
	double vy;
	double vz;

	int sector_e;
	double p_e; 		
	double theta_e; 	
	double phi_e;		
	double p_n;		
	double theta_n;	
	double phi_n;		
	double E_n;
	double Q2;		
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
	double theta_nq;
	double phi_nq;
	double CosTheta_nq;
	double Xp;
	double Wp;
	double As;


};

#endif

