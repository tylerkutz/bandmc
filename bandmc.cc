#include "ElectronFiducial.hh"
#include "Smearing.hh"
#include "GenTree.hh"
#include "IO.hh"
#include "BAND.hh"
#include "RadGen.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

	// Read in arguments
	if (argc != 2)
	{
		cerr << "Wrong number of arguments. Instead use: [/path/to/input/file]\n";
		exit(-1);
	}

	double Me = 0.511e-3;
	double Mp = 0.938272;
	double eBeam = 10.6;
	float ebf = 10.6;

	TFile* inputFile  = new TFile(argv[1]);
	TTree* mcTree=(TTree*)inputFile->Get("MCout");

	Gen_Event* event = new Gen_Event();
	Gen_Particle* electron = new Gen_Particle(); 
	Gen_Particle* neutron = new Gen_Particle(); 

	ElectronFiducial* fFiducial = new ElectronFiducial("upper_momentum_fit.dat", "lower_momentum_fit.dat");
	BAND* fBAND = new BAND();
	RadGen* fRad = new RadGen();
	Smearing* fSmear = new Smearing();
	IO* fIO = new IO("bandmc_out.root");

	mcTree->SetBranchAddress("event", &event);

//	TH2F* hist = new TH2F("hist","", 420, -220, 200, 160, 0, 40);

	int nEvents = mcTree->GetEntries();
	nEvents = 1000;

	for (int i = 0; i<nEvents; i++) {

		if(i%100000==0) {
			cout << i << "/" << nEvents << endl;
		}

		mcTree->GetEntry(i);

		// Electron info
		electron = &event->particles[0];
		double px_e = electron->momentum.x();
		double py_e = electron->momentum.y();
		double pz_e = electron->momentum.z();
		double p_e = electron->momentum.Mag(); 
		double theta_e = electron->momentum.Theta()*TMath::RadToDeg();
		double phi_e = electron->momentum.Phi()*TMath::RadToDeg();
		double E_e = sqrt(p_e*p_e + Me*Me );

		// Neutron info
		neutron = &event->particles[1];
		double p_n = neutron->momentum.Mag();
		double theta_n = neutron->momentum.Theta()*TMath::RadToDeg();
		double phi_n = neutron->momentum.Phi()*TMath::RadToDeg();

		// Event info
		TVector3 E0 = TVector3(0.,0.,eBeam);
		TVector3 q = electron->momentum - E0;
		double Q2 = q.Mag()*q.Mag();
		double nu = eBeam - E_e;
		double xB = Q2/(2*Mp*nu);
		
		// Radiate
		float vpgen[4] = {(float)q.x(), (float)q.y(), (float)q.z(), (float)nu};
		float vprad[4];
		float rprad[4];
		float Q2true;
		float Utrue;
		float weight;
		fRad->Radiate(&ebf, vpgen, vprad, rprad, &Q2true, &Utrue, &weight); 

		// Smear detected electron
		TLorentzVector e4;
		e4.SetPxPyPzE(px_e, py_e, pz_e, E_e);
		fSmear->smear(&e4, 0);

		// Check acceptance for electron, neutron
		double acc_e = fFiducial->GetElectronAcceptance(theta_e, phi_e, p_e);
		double acc_n = fBAND->GetNeutronAcceptance(theta_n, phi_n, 0.);

		if (acc_e > 0.5 && acc_n > 0.5) {

			fIO->p_e = e4.P();
			fIO->theta_e = e4.Theta()*TMath::RadToDeg();
			fIO->phi_e =  e4.Phi()*TMath::RadToDeg();
			fIO->p_n = p_n;
			fIO->theta_n = theta_n;
			fIO->phi_n = phi_n;
	
			fIO->FillTree();

		}


	}

	fIO->WriteTree();

	return 1;
}


