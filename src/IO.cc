#include "IO.hh"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <unistd.h>

using namespace std;

IO::IO(TString filename) {

	fFile = new TFile(filename, "RECREATE");
	
	InitializeTree();

	nMult = 1;

}

IO::~IO() {


}

void IO::InitializeTree() {

	fTree = new TTree("T", "BAND fast Monte Carlo");	
	
	fTree->SetMaxTreeSize(1900000000);  // 1.9 GB

	fTree->Branch("p_e",		&p_e, 		"p_e/D");	
	fTree->Branch("theta_e",	&theta_e, 	"theta_e/D");	
	fTree->Branch("phi_e",		&phi_e,		"phi_e/D");
	fTree->Branch("p_n",		&p_n,		"p_n/D");
	fTree->Branch("theta_n",	&theta_n,	"theta_n/D");
	fTree->Branch("phi_n",		&phi_n,		"phi_n/D");
	fTree->Branch("Q2",		&Q2,		"Q2/D");
	fTree->Branch("xB",		&xB,		"xB/D");		
	fTree->Branch("W2",		&W2,		"W2/D");
	fTree->Branch("Ebeam",		&Ebeam,		"Ebeam/D");
	fTree->Branch("gated_charge",	&gated_charge,	"gated_charge/D");
	fTree->Branch("livetime",	&livetime,	"livetime/D");
	fTree->Branch("starttime",	&starttime,	"starttime/D");
	fTree->Branch("ePid",		&ePid,		"ePid/I");
	fTree->Branch("eCharge",	&eCharge,	"eCharge/D");
	fTree->Branch("eStatus",	&eStatus,	"eStatus/I");
	fTree->Branch("eTime",		&eTime,		"eTime/D");
	fTree->Branch("eBeta",		&eBeta,		"eBeta/D");
	fTree->Branch("eChi2pid",	&eChi2pid,	"eChi2pid");
	fTree->Branch("E_tot",		&E_tot,		"E_tot/D");
	fTree->Branch("E_pcal",		&E_pcal,	"E_pcal/D");
	fTree->Branch("t_e",		&t_e,		"t_e/D");
	fTree->Branch("dL_e",		&dL_e,		"dL_e/D");
	fTree->Branch("lU",		&lU,		"lU/D");
	fTree->Branch("lV",		&lV,		"lV/D");
	fTree->Branch("lW",		&lW,		"lW/D");
	fTree->Branch("e_vtx",		&e_vtx,		"e_vtx/D");
	fTree->Branch("e_vty",		&e_vty,		"e_vty/D");
	fTree->Branch("e_vtz",		&e_vtz,		"e_vtz/D");
	fTree->Branch("q",		&q,		"q/D");
	fTree->Branch("theta_q",	&theta_q,	"theta_q/D");
	fTree->Branch("phi_q",		&phi_q,		"phi_q/D");
	fTree->Branch("nu",		&nu,		"nu/D");
	fTree->Branch("nMult",		&nMult,		"nMult/I");
	fTree->Branch("barID",		&barID,		"barID/D");
	fTree->Branch("dL_n",		&dL_n,		"dL_n/D");
	fTree->Branch("nTime",		&nTime,		"nTime/D");
	fTree->Branch("nEdep",		&nEdep,		"nEdep/D");
	
}


void IO::FillTree() {

	fTree->Fill();

}

void IO::WriteTree() {

	fFile->cd();
	fTree->Write("T", TObject::kOverwrite);

	fTree->ResetBranchAddresses();
	delete fTree;
	fTree = NULL;

	fFile->Close();

	delete fFile;
	fFile = NULL;


}

void IO::ResetBranches() {

	p_e = 0; 		
	theta_e = 0; 	
	phi_e = 0;		
	p_n = 0;		
	theta_n = 0;	
	phi_n = 0;		
	Q2 = 0;		
	xB = 0;		
	W2 = 0;		
	Ebeam = 0;		
	gated_charge = 0;	
	livetime = 0;	
	starttime = 0;	
	ePid = 0;		
	eCharge = 0;	
	eStatus = 0;	
	eTime = 0;		
	eBeta = 0;		
	eChi2pid = 0;	
	E_tot = 0;		
	E_pcal = 0;	
	t_e = 0;		
	dL_e = 0;		
	lU = 0;		
	lV = 0;		
	lW = 0;		
	e_vtx = 0;		
	e_vty = 0;		
	e_vtz = 0;		
	q = 0;		
	theta_q = 0;	
	phi_q = 0;		
	nu = 0;		
	nMult = 0;		
	barID = 0;		
	dL_n = 0;		
	nTime = 0;		
	nEdep = 0;		

}
