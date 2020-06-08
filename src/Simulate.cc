#include "Simulate.hh"
#include "RadGen.hh"
#include "Smearing.hh"
#include "IO.hh"
#include "ElectronFiducial.hh"
#include "BAND.hh"

#include "TVector3.h"
#include "TRandom3.h"

#include <iostream>
using namespace std;

Simulate::Simulate(int do_radiate, int do_smear) {

	doRadiation = do_radiate;
	doSmearing = do_smear;

	fRad = new RadGen();
	fSmear = new Smearing();
	fFiducial = new ElectronFiducial("./dat/upper_momentum_fit.dat", "./dat/lower_momentum_fit.dat");
	fBAND = new BAND();
	
	fRand = new TRandom3(0);

	Me = 0.511e-3;
	Mp = 0.938272;
	Mn = 0.939565;
	eBeam = 10.6;
	c = 2.9979e8;
	MD = 1.8756;

	nSim = 0;
	nBG = 0;

	targetL = 5.;

	barThickness = 7.2;
	// Distance from ideal target center to upstream face of BAND 
	zUp = ((302.69+302.69+302.57+302.64)/4.) ; // [cm]
	// Distance from ideal target center to downstream face of BAND
	zDown = zUp - 5.*barThickness;


}

Simulate::~Simulate() {

	delete fRad;
	delete fSmear;
	delete fFiducial;
	delete fBAND;
	delete fIO;
	delete fRand;

}

int Simulate::SimulateEvent(TVector3 electron, TVector3 neutron) {

		fIO->ResetBranches();

		bg = 0;

		// Vertex info
		double vz = fRand->Uniform(-targetL/2., targetL/2.);

		// Electron info
		double px_e = electron.x();
		double py_e = electron.y();
		double pz_e = electron.z();
		p_e = electron.Mag(); 
		theta_e = electron.Theta()*TMath::RadToDeg();
		phi_e = electron.Phi()*TMath::RadToDeg();
		double E_e = sqrt(p_e*p_e + Me*Me );
		eTime = 0;

		// Neutron info
		p_n = neutron.Mag();
		theta_n = neutron.Theta()*TMath::RadToDeg();
		phi_n = neutron.Phi()*TMath::RadToDeg();
		double E_n = sqrt(Mn*Mn + p_n*p_n );

		double nGamma2 = 1. + (p_n/Mn)*(p_n/Mn);
		nBeta = sqrt(1. - 1./nGamma2);
		double nv = nBeta*c;
		double dZ_n = -fRand->Uniform(zDown + vz, zUp + vz);  
		dL_n = dZ_n/cos(neutron.Theta()); 	
		nTime = 1.e9*(dL_n/100.)/nv;		// convert dL_n to meters first, give time in ns 

		// Unradiated photon info
		TVector3 k0 = TVector3(0.,0.,sqrt(eBeam*eBeam - Me*Me));
		TVector3 q = k0 - electron;
		nu = eBeam - E_e;	
		Q2 = q.Mag()*q.Mag() - nu*nu;
		theta_q = q.Theta();
		phi_q = q.Phi();

		// nq vectors
		TVector3 norm_scatter = q.Cross( k0 );
		norm_scatter = norm_scatter.Unit();
		TVector3 norm_reaction = q.Cross( neutron );
		norm_reaction   = norm_reaction.Unit();

		double phi_nq   = norm_scatter.Angle( norm_reaction );
		double theta_nq = neutron.Angle( q );
	
		CosTheta_nq = cos(theta_nq);
		TVector3 direction = norm_scatter.Cross(norm_reaction);
	
		if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
		}
		else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
			phi_nq *= (-1);
		}
		double W_primeSq = MD*MD - Q2 + Mn*Mn + 2.*MD*(nu-E_n) - 2.*nu*E_n + 2.*q.Mag()*p_n*cos(theta_nq);
		Wp = sqrt(W_primeSq);
		Xp = Q2/(2.*( nu*(MD-E_n) + p_n*q.Mag()*CosTheta_nq));
		As = (E_n - p_n*CosTheta_nq)/Mn;

		// Radiate photon
		float vpgen[4] = {(float)q.x(), (float)q.y(), (float)q.z(), (float)nu};
		float vprad[4];
		float rprad[4];
		float Q2true;
		float Utrue;
		float weight;

		if (doRadiation) {
			float ebf = (float)eBeam;	
			fRad->Radiate(&ebf, vpgen, vprad, rprad, &Q2true, &Utrue, &weight); 
			Q2 = Q2true;
			q = TVector3(vprad[0], vprad[1], vprad[2]);
			theta_q = q.Theta();
			phi_q = q.Phi();			
		}

		// Smear detected electron
//		cout << "before " << p_e << endl;
		e4.SetPxPyPzE(px_e, py_e, pz_e, E_e);
		if(doSmearing) {
			fSmear->smear(&e4, 0);
		}
		p_e = e4.Vect().Mag();
//		cout << "after " << p_e << endl;
		E_e = e4.E();
		theta_e = e4.Vect().Theta()*TMath::RadToDeg();
		phi_e = e4.Vect().Phi()*TMath::RadToDeg();
		nu = eBeam - E_e;

		// After radiation and smearing, calculate remaining kinematic variables
		xB = Q2/(2*Mp*nu);
		W2 =  Mp*Mp + 2.*nu*Mp - Q2;
		double p_mc2 = (p_e/Me)*(p_e/Me);
	      	eBeta = sqrt(p_mc2 / (1. + p_mc2));	


		// Check acceptance for electron, neutron
		double acc_e = fFiducial->GetElectronAcceptance(theta_e, phi_e, p_e);
		double acc_n = fBAND->GetNeutronAcceptance(theta_n, phi_n, p_n, vz);

		if (acc_e > 0.5 && acc_n > 0.5) {
			
			SetData();
			fIO->FillTree();
			nSim++;

		}

	return nSim;

}

int Simulate::SimulateBackground(TVector3 electron, TVector3 neutron) {

		fIO->ResetBranches();

		bg = 1;

		// Vertex info
		double vz = fRand->Uniform(-targetL/2., targetL/2.);

		// Electron info
		double px_e = electron.x();
		double py_e = electron.y();
		double pz_e = electron.z();
		p_e = electron.Mag(); 
		theta_e = electron.Theta()*TMath::RadToDeg();
		phi_e = electron.Phi()*TMath::RadToDeg();
		double E_e = sqrt(p_e*p_e + Me*Me );
		eTime = 0;

		// Neutron info
		p_n = neutron.Mag();
		theta_n = neutron.Theta()*TMath::RadToDeg();
		phi_n = neutron.Phi()*TMath::RadToDeg();
		double E_n = sqrt(Mn*Mn + p_n*p_n );

		double nGamma2 = 1. + (p_n/Mn)*(p_n/Mn);
		nBeta = sqrt(1. - 1./nGamma2);
		double nv = nBeta*c;
		double dZ_n = -fRand->Uniform(zDown + vz, zUp + vz);  
		dL_n = dZ_n/cos(neutron.Theta()); 	
		nTime = fRand->Gaus(fRand->Uniform(-20., 340.), 5.);

		// Unradiated photon info
		TVector3 k0 = TVector3(0.,0.,sqrt(eBeam*eBeam - Me*Me));
		TVector3 q = k0 - electron;
		nu = eBeam - E_e;	
		Q2 = q.Mag()*q.Mag() - nu*nu;
		theta_q = q.Theta();
		phi_q = q.Phi();

		// nq vectors
		TVector3 norm_scatter = q.Cross( k0 );
		norm_scatter = norm_scatter.Unit();
		TVector3 norm_reaction = q.Cross( neutron );
		norm_reaction   = norm_reaction.Unit();

		double phi_nq   = norm_scatter.Angle( norm_reaction );
		double theta_nq = neutron.Angle( q );
	
		CosTheta_nq = cos(theta_nq);
		TVector3 direction = norm_scatter.Cross(norm_reaction);
	
		if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
		}
		else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
			phi_nq *= (-1);
		}
		double W_primeSq = MD*MD - Q2 + Mn*Mn + 2.*MD*(nu-E_n) - 2.*nu*E_n + 2.*q.Mag()*p_n*cos(theta_nq);
		Wp = sqrt(W_primeSq);
		Xp = Q2/(2.*( nu*(MD-E_n) + p_n*q.Mag()*CosTheta_nq));
		As = (E_n - p_n*CosTheta_nq)/Mn;

		// Radiate photon
		float vpgen[4] = {(float)q.x(), (float)q.y(), (float)q.z(), (float)nu};
		float vprad[4];
		float rprad[4];
		float Q2true;
		float Utrue;
		float weight;

		if (doRadiation) {
			float ebf = (float)eBeam;	
			fRad->Radiate(&ebf, vpgen, vprad, rprad, &Q2true, &Utrue, &weight); 
			Q2 = Q2true;
			q = TVector3(vprad[0], vprad[1], vprad[2]);
			theta_q = q.Theta();
			phi_q = q.Phi();			
		}

		// Smear detected electron
//		cout << "before " << p_e << endl;
		e4.SetPxPyPzE(px_e, py_e, pz_e, E_e);
		if(doSmearing) {
			fSmear->smear(&e4, 0);
		}
		p_e = e4.Vect().Mag();
//		cout << "after " << p_e << endl;
		E_e = e4.E();
		theta_e = e4.Vect().Theta()*TMath::RadToDeg();
		phi_e = e4.Vect().Phi()*TMath::RadToDeg();
		nu = eBeam - E_e;

		// After radiation and smearing, calculate remaining kinematic variables
		xB = Q2/(2*Mp*nu);
		W2 =  Mp*Mp + 2.*nu*Mp - Q2;
		double p_mc2 = (p_e/Me)*(p_e/Me);
	      	eBeta = sqrt(p_mc2 / (1. + p_mc2));	


		// Check acceptance for electron, neutron
		double acc_e = fFiducial->GetElectronAcceptance(theta_e, phi_e, p_e);
		double acc_n = fBAND->GetNeutronAcceptance(theta_n, phi_n, p_n, vz);

		if (acc_e > 0.5 && acc_n > 0.5) {
			
			SetData();
			fIO->FillTree();
			nBG++;

		}

	return nBG;

}

void Simulate::SetData() {

	fIO->p_e = e4.P();
	fIO->theta_e = e4.Theta()*TMath::RadToDeg();
	fIO->phi_e =  e4.Phi()*TMath::RadToDeg();
	fIO->p_n = p_n;
	fIO->theta_n = theta_n;
	fIO->phi_n = phi_n;
	fIO->Q2 = Q2;
	fIO->xB = xB;
	fIO->W2 = W2;
	fIO->eBeta = eBeta;
	fIO->E_tot = E_e;
	fIO->theta_q = theta_q;
	fIO->phi_q = phi_q;
	fIO->nTime = nTime;	
	fIO->dL_n = dL_n;
	fIO->nBeta = nBeta;
	fIO->bg = bg;
	fIO->CosTheta_nq = CosTheta_nq;
	fIO->Xp = Xp;
	fIO->Wp = Wp;
	fIO->As = As;

}

