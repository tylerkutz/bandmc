#include "IO.hh"

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include <unistd.h>

using namespace std;

IO::IO(TString filename) {

	fFile = new TFile(filename, "RECREATE");

	// These things will not change for each MC run or event
	fRunNum = 1;
	fEbeam = 10.6; // GeV
	fGatedCharge = 1.;
	fLivetime = 1.;
	fStartTime = 0.;
	fCurrent = 1.;
	fnMult = 1;

	InitializeTree();

}

IO::~IO() {


}

void IO::InitializeTree() {

	fTree = new TTree("T", "BAND fast Monte Carlo");	
	
	fTree->SetMaxTreeSize(1900000000);  // 1.9 GB

        fTree->Branch("Runno",         	&fRunNum);
        fTree->Branch("Ebeam",         	&fEbeam);
        fTree->Branch("gated_charge",  	&fGatedCharge);
        fTree->Branch("livetime",      	&fLivetime);
        fTree->Branch("starttime",     	&fStartTime);
        fTree->Branch("current",       	&fCurrent);
        fTree->Branch("nMult",         	&fnMult);
        fTree->Branch("nHits",          &fBANDHit);
        fTree->Branch("eHit",         	&fCLASHit);
        fTree->Branch("tag",         	&fTagHit);
	fTree->Branch("bg",		&bg);
	fTree->Branch("radweight",	&radweight);

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

void IO::ClearEvent() {

	fCLASHit.Clear();
	fBANDHit->Clear();
	fTagHit->Clear();
	bg = -1;
	radweight = 1.;

}
