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

	int numargs = 3;
	// Read in arguments
	if (argc < numargs) {
		cerr << "Wrong number of arguments. Instead use: [input file] [output file]\n";
		exit(-1);
	}

	bool doSmearing = false;
       	bool isInclusive = false;	
	bool doBackground = false;
	int nBackground = 0;

	int c;
	while ((c = getopt (argc-numargs+1, &argv[numargs-1], "B:SI")) != -1) { 

		switch(c) {

			case 'B':
				doBackground = true;
				nBackground = atoi(optarg);
				break;
			case 'S':
				doSmearing = true;
				break;
			case 'I': 
				isInclusive = true;
				break;
			default:
				break;
				//do nothing
		}
	}		
	

	Gen_Event* event = new Gen_Event();
	Gen_Event* backg = new Gen_Event();


//	DISGenerator* fDIS = new DISGenerator(doRadiation); 
	TFile* genFile = new TFile(argv[1]);
	TTree* genTree = (TTree*)genFile->Get("T");
	double pe[3], pn[3];
	double radweight;
	genTree->SetBranchAddress("pe", pe);
//	genTree->SetBranchAddress("radweight", &radweight);
	if(!isInclusive) {
		genTree->SetBranchAddress("pn", pn);
	}
	int nEvents = genTree->GetEntries();

	cout << "Smearing:  " << doSmearing << endl;
	cout << "Inclusive: " << isInclusive << endl;

	RandomGenerator* fBG = new RandomGenerator();
	Simulate* fSim = new Simulate(doSmearing, isInclusive);

	IO* fIO = new IO(argv[2]);

	fSim->SetIO(fIO);

	int nSim = 0;
	int nBG = 0;

	for(int i = 0; i < nEvents; i++) {

		if( i%100000 == 0 ) {
			cout << "Processed " << i << "/" << nEvents << endl;
		}

		genTree->GetEntry(i);
		nSim = fSim->SimulateEvent(TVector3(pe[0], pe[1], pe[2]), TVector3(pn[0], pn[1], pn[2]), radweight);
		

	}


	while (nBG < nBackground) {

		if(nBG%100000==0) {
			cout << "Background: " << nBG << "/" << nBackground << endl;
		}
		fBG->GenerateRandom(backg);
		nBG = fSim->SimulateBackground(backg->particles[0].momentum, backg->particles[1].momentum);

	}




	fIO->WriteTree();

	return 1;
}


