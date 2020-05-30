#ifndef IO_HH
#define IO_HH

class TString;
class TFile;
class TTree;
class TVector3;

class IO{

public:
	IO(TString);
	~IO();

	void FillTree();
	void WriteTree();
	void ResetBranches();


	double p_e; 		
	double theta_e; 	
	double phi_e;		
	double p_n;		
	double theta_n;	
	double phi_n;		
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
	double nEdep;		


private:

	TFile* fFile;
	TTree* fTree;

	void InitializeTree();

};

#endif

