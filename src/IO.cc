#include "IO.hh"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <unistd.h>

using namespace std;

IO::IO(TString filename) {

	fFile = new TFile(filename, "RECREATE");
	
	InitializeTree();

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

