#include "Radiation.hh"

#include "TGraph2D.h"
#include "TFile.h"

#include <fstream>
#include <vector>

using namespace std;

Radiation::Radiation() {

	LoadGraph("dat/RC_graph.root");

}

Radiation::~Radiation() {

}


double Radiation::GetRadiativeCorrection(double x, double Q2) {

	double RC = fRadGraph->Interpolate(x, Q2);

	if(RC < 0.001) {
		RC = 1.0;
	}

	return RC;

}

void Radiation::LoadGraph(TString radfilename) {

	TFile* radFile = new TFile(radfilename);

	fRadGraph = (TGraph2D*)radFile->Get("rcgr");

}



