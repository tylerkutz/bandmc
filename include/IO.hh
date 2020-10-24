#ifndef IO_HH
#define IO_HH

#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "TClonesArray.h"

class TFile;
class TTree;

class IO{

public:
	IO(TString);
	~IO();

	void FillTree();
	void WriteTree();
	void ClearEvent();

	int fRunNum;              
        double fEbeam;           
        double fGatedCharge;    
        double fLivetime;        
        double fStartTime;      
        double fCurrent;     
        //      Neutron info:
        int fnMult;              
//	bandhit fBANDHit;
//	taghit fTagHit;
	TClonesArray* fBANDHit = new TClonesArray("bandhit");
        TClonesArray &saveHit = *fBANDHit;
	TClonesArray* fTagHit = new TClonesArray("taghit");
        TClonesArray &saveTag = *fTagHit;
	clashit fCLASHit;
	double radweight;
	int bg;

private:

	TFile* fFile;
	TTree* fTree;

	void InitializeTree();

};

#endif

