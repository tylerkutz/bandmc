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

using namespace std;

int main(int argc, char *argv[]) {

	// Read in arguments
	if (argc != 5)
	{
		cerr << "Wrong number of arguments. Instead use: [N events] [N background] [radiation on/off {1,0}] [smearing on/off {1,0}]\n";
		exit(-1);
	}


	int nEvents = atoi(argv[1]);
	int nBackg = atoi(argv[2]);

	int doRadiation = atoi(argv[3]);
	int doSmearing = atoi(argv[4]);

	Gen_Event* event = new Gen_Event();
	Gen_Event* backg = new Gen_Event();

	IO* fIO = new IO(Form("bandmc_out_rad%i_smear%i.root", doRadiation, doSmearing));

	DISGenerator* fDIS = new DISGenerator(doRadiation); 

	RandomGenerator* fBG = new RandomGenerator();
	Simulate* fSim = new Simulate(doSmearing);

	fSim->SetIO(fIO);

	int nSim = 0;
	int nBG = 0;

	while (nSim < nEvents) {

		if(nSim%100==0) {
			cout << "Signal: " << nSim << "/" << nEvents << endl;
		}
		fDIS->GenerateEvent(event);		
		nSim = fSim->SimulateEvent(event->particles[0].momentum, event->particles[1].momentum);

	}


	while (nBG < nBackg) {

		if(nBG%100==0) {
			cout << "Background: " << nBG << "/" << nBackg << endl;
		}
		fBG->GenerateRandom(backg);
		nBG = fSim->SimulateBackground(backg->particles[0].momentum, backg->particles[1].momentum);

	}




	fIO->WriteTree();

	return 1;
}


