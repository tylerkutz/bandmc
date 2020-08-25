#include "DISGenerator.hh"
#include "RandomGenerator.hh"
#include "GenTree.hh"
#include "IO.hh"
#include "Simulate.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TMath.h"

#include <iostream>
#include <unistd.h>

using namespace std;

int main(int argc, char *argv[]) {

	int numargs = 2;

	// Read in arguments
	if (argc < numargs) {
		cerr << "Wrong number of arguments. Instead use: [input file]\n";
		exit(-1);
	}

	bool doBackground = false;
	bool doSmearing = false;
       	bool isInclusive = false;	

	int nBackground = 0;

	int c;
	while ((c = getopt (argc-numargs+1, &argv[numargs-1], "b:si")) != -1) {
		
		switch(c) {

			case 'b':
				doBackground = true;
				nBackground = atoi(optarg);
				break;
			case 's':
				doSmearing = true;
				break;
			case 'i': 
				isInclusive = true;
			default:
				abort();

		}
	}

	Gen_Event* event = new Gen_Event();
	Gen_Event* backg = new Gen_Event();

	IO* fIO = new IO("bandmc_out.root");

//	DISGenerator* fDIS = new DISGenerator(doRadiation); 
	TFile* genFile = new TFile(argv[1]);
	TTree* genTree = (TTree*)genFile->Get("T");
	double pe[3], pn[3];
	genTree->SetBranchAddress("pe", &pe);
	genTree->SetBranchAddress("pn", &pn);
	int nEvents = genTree->GetEntries();

	RandomGenerator* fBG = new RandomGenerator();
	Simulate* fSim = new Simulate(doSmearing, isInclusive);

	fSim->SetIO(fIO);

	int nSim = 0;
	int nBG = 0;

	for(int i = 0; i < nEvents; i++) {

//		if(nSim%100==0) {
//			cout << "Signal: " << nSim << "/" << nEvents << endl;
//		}
//		fDIS->GenerateEvent(event);		
//		nSim = fSim->SimulateEvent(event->particles[0].momentum, event->particles[1].momentum);

		genTree->GetEntry(i);
		nSim = fSim->SimulateEvent(TVector3(pe[0], pe[1], pe[2]), TVector3(pn[0], pn[1], pn[2]));
		

	}


	while (nBG < nBackground) {

		if(nBG%100==0) {
			cout << "Background: " << nBG << "/" << nBackground << endl;
		}
		fBG->GenerateRandom(backg);
		nBG = fSim->SimulateBackground(backg->particles[0].momentum, backg->particles[1].momentum);

	}




	fIO->WriteTree();

	return 1;
}


