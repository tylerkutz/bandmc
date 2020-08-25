#include "DISGenerator.hh"
#include "GenTree.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

	// Read in arguments
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead use: [N events] [radiation on/off {1,0}] \n";
		exit(-1);
	}


	int nEvents = atoi(argv[1]);
	int doRadiation = atoi(argv[2]);

	Gen_Event* event = new Gen_Event();

	DISGenerator* fDIS = new DISGenerator(doRadiation); 

	double pe[3];
	double pn[3];

	TFile* fFile = new TFile("band_wim_tagged.root", "RECREATE");	
	TTree* fTree = new TTree("T", "BAND fast Monte Carlo");

	fTree->Branch("pe", pe, "pe[3]/D");
	fTree->Branch("pn", pn, "pn[3]/D");

	for(int i = 0; i < nEvents; i++) {

		fDIS->GenerateEvent(event);

		pe[0] = event->particles[0].momentum.x();
		pe[1] = event->particles[0].momentum.y();
		pe[2] = event->particles[0].momentum.z();

		pn[0] = event->particles[1].momentum.x();
		pn[1] = event->particles[1].momentum.y();
		pn[2] = event->particles[1].momentum.z();
	
		fTree->Fill();

	}

        fFile->cd();
        fTree->Write("T", TObject::kOverwrite);

        fFile->Close();

        delete fFile;
        fFile = NULL;
	
	return 1;
}


