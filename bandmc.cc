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
#include "TMath.h"

#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

	// Read in arguments
	if (argc != 4)
	{
		cerr << "Wrong number of arguments. Instead use: [/path/to/input/file] [radiation on/off {1,0}] [smearing on/off {1,0}]\n";
		exit(-1);
	}

	double Me = 0.511e-3;
	double Mp = 0.938272;
	double eBeam = 10.6;

	TFile* inputFile  = new TFile(argv[1]);
	TTree* mcTree=(TTree*)inputFile->Get("MCout");

	int do_radiate = atoi(argv[2]);
	int do_smearing = atoi(argv[3]);

	Gen_Event* event = new Gen_Event();
	Gen_Particle* electron = new Gen_Particle(); 
	Gen_Particle* neutron = new Gen_Particle(); 

	ElectronFiducial* fFiducial = new ElectronFiducial("upper_momentum_fit.dat", "lower_momentum_fit.dat");
	BAND* fBAND = new BAND();
	RadGen* fRad = new RadGen();
	Smearing* fSmear = new Smearing();
	IO* fIO = new IO(Form("bandmc_out_rad%i_smear%i.root", do_radiate, do_smearing));

	mcTree->SetBranchAddress("event", &event);

//	TH2F* hist = new TH2F("hist","", 420, -220, 200, 160, 0, 40);

	int nEvents = mcTree->GetEntries();
	nEvents = 1000000;

	for (int i = 0; i<nEvents; i++) {

		if(i%100==0) {
			cout << i << "/" << nEvents << endl;
		}

		fIO->ResetBranches();

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

		// Unradiated photon info
		TVector3 k0 = TVector3(0.,0.,sqrt(eBeam*eBeam - Me*Me));
		TVector3 q = electron->momentum - k0;
		double nu = eBeam - E_e;	
		double Q2 = q.Mag()*q.Mag() - nu*nu;
		double theta_q = q.Theta();
		double phi_q = q.Phi();

		// Radiate photon
		float vpgen[4] = {(float)q.x(), (float)q.y(), (float)q.z(), (float)nu};
		float vprad[4];
		float rprad[4];
		float Q2true;
		float Utrue;
		float weight;

		if (do_radiate) {
			float ebf = (float)eBeam;	
			fRad->Radiate(&ebf, vpgen, vprad, rprad, &Q2true, &Utrue, &weight); 
			Q2 = Q2true;
			q = TVector3(vprad[0], vprad[1], vprad[2]);
			theta_q = q.Theta();
			phi_q = q.Phi();			
		}

		// Smear detected electron
//		cout << "before " << p_e << endl;
		TLorentzVector e4;
		e4.SetPxPyPzE(px_e, py_e, pz_e, E_e);
		if(do_smearing) {
			fSmear->smear(&e4, 0);
		}
		p_e = e4.Vect().Mag();
//		cout << "after " << p_e << endl;
		E_e = e4.E();
		theta_e = e4.Vect().Theta()*TMath::RadToDeg();
		phi_e = e4.Vect().Phi()*TMath::RadToDeg();
		nu = eBeam - E_e;

		// After radiation and smearing, calculate remaining kinematic variables
		double xB = Q2/(2*Mp*nu);
		double W2 =  Mp*Mp + 2.*nu*Mp - Q2;
		double p_mc2 = (p_e/Me)*(p_e/Me);
	      	double beta = sqrt(p_mc2 / (1. + p_mc2));	

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
			fIO->Q2 = Q2;
			fIO->xB = xB;
			fIO->W2 = W2;
			fIO->eBeta = beta;
			fIO->E_tot = E_e;
			fIO->theta_q = theta_q;
			fIO->phi_q = phi_q;
			fIO->FillTree();

		}


	}

	fIO->WriteTree();

	return 1;
}


