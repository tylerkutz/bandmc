#include "Simulate.hh"
#include "IO.hh"
#include "CLAS.hh"
#include "BAND.hh"

#include "TVector3.h"
#include "TRandom3.h"


#include <iostream>
using namespace std;

Simulate::Simulate(int do_smear) {

	doSmearing = do_smear;

	fCLAS = new CLAS();
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
	targetR = 1.;

	barThickness = 7.2;
	// Distance from ideal target center to upstream face of BAND 
	zUp = ((302.69+302.69+302.57+302.64)/4.) ; // [cm]
	// Distance from ideal target center to downstream face of BAND
	zDown = zUp - 5.*barThickness;


}

Simulate::~Simulate() {

	delete fCLAS;
	delete fBAND;
	delete fIO;
	delete fRand;

}

int Simulate::SimulateEvent(TVector3 electron, TVector3 neutron) {

		fIO->ClearEvent();
		
		bg = 0;

		// Vertex
		vx = fRand->Uniform(targetR) * cos(fRand->Uniform(360.*TMath::DegToRad()));
		vy = fRand->Uniform(targetR) * sin(fRand->Uniform(360.*TMath::DegToRad()));
		vz = fRand->Uniform(-targetL/2., targetL/2.);

		// Create Lorentz vector required for electron smearing
		double px_e = electron.x();
		double py_e = electron.y();
		double pz_e = electron.z();
		p_e = electron.Mag(); 
		E_e = sqrt(p_e*p_e + Me*Me);
		e4.SetPxPyPzE(px_e, py_e, pz_e, E_e);

		// Electron and neutron smearing
		if(doSmearing) {
			fCLAS->Smear(&e4, 0);
			fBAND->Smear(&neutron);
		}

		// Electron variables
		p_e = e4.Vect().Mag();
		E_e = e4.E();
		theta_e = e4.Vect().Theta()*TMath::RadToDeg();
		phi_e = e4.Vect().Phi()*TMath::RadToDeg();
		double p_mc2 = (p_e/Me)*(p_e/Me);
	      	eBeta = sqrt(p_mc2 / (1. + p_mc2));	
		eTime = 0;

		// Neutron variables
		p_n = neutron.Mag();
		theta_n = neutron.Theta()*TMath::RadToDeg();
		phi_n = neutron.Phi()*TMath::RadToDeg();
		E_n = sqrt(Mn*Mn + p_n*p_n );
		double nGamma2 = 1. + (p_n/Mn)*(p_n/Mn);
		nBeta = sqrt(1. - 1./nGamma2);
		double nv = nBeta*c;
		double dZ_n = -fRand->Uniform(zDown + vz, zUp + vz);  
		dL_n = dZ_n/cos(neutron.Theta()); 	
		nTime = 1.e9*(dL_n/100.)/nv;		// convert dL_n to meters first, give time in ns 
		nEdep = fBAND->GetEdep(p_n);

		// Photon variables
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
		phi_nq   = norm_scatter.Angle( norm_reaction );
		theta_nq = neutron.Angle( q );
		CosTheta_nq = cos(theta_nq);
		TVector3 direction = norm_scatter.Cross(norm_reaction);
	
		if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
		}
		else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
			phi_nq *= (-1);
		}

		// Remaining kinematic variables
		xB = Q2/(2*Mp*nu);
		W2 =  Mp*Mp + 2.*nu*Mp - Q2;
		double W_primeSq = MD*MD - Q2 + Mn*Mn + 2.*MD*(nu-E_n) - 2.*nu*E_n + 2.*q.Mag()*p_n*cos(theta_nq);
		Wp = sqrt(W_primeSq);
		Xp = Q2/(2.*( nu*(MD-E_n) + p_n*q.Mag()*CosTheta_nq));
		As = (E_n - p_n*CosTheta_nq)/Mn;

		// Check acceptance for electron, neutron
		sector_e = fCLAS->GetElectronAcceptance(theta_e, phi_e, p_e);
		double acc_n = fBAND->GetNeutronAcceptance(theta_n, phi_n, p_n, vz);

		if (sector_e > 0 && acc_n > 0.5) {


			SetEventData();
			fIO->FillTree();
			nSim++;

		}

	return nSim;

}

int Simulate::SimulateBackground(TVector3 electron, TVector3 neutron) {

		fIO->ClearEvent();

		bg = 1;

		// Vertex
		vx = fRand->Uniform(targetR) * cos(fRand->Uniform(360.*TMath::DegToRad()));
		vy = fRand->Uniform(targetR) * sin(fRand->Uniform(360.*TMath::DegToRad()));
		vz = fRand->Uniform(-targetL/2., targetL/2.);

		// Create Lorentz vector required for electron smearing
		double px_e = electron.x();
		double py_e = electron.y();
		double pz_e = electron.z();
		p_e = electron.Mag(); 
		E_e = sqrt(p_e*p_e + Me*Me);
		e4.SetPxPyPzE(px_e, py_e, pz_e, E_e);

		// Electron and neutron smearing
		if(doSmearing) {
			fCLAS->Smear(&e4, 0);
			fBAND->Smear(&neutron);
		}

		// Electron variables
		p_e = e4.Vect().Mag();
		E_e = e4.E();
		theta_e = e4.Vect().Theta()*TMath::RadToDeg();
		phi_e = e4.Vect().Phi()*TMath::RadToDeg();
		double p_mc2 = (p_e/Me)*(p_e/Me);
	      	eBeta = sqrt(p_mc2 / (1. + p_mc2));	
		eTime = 0;

		// Use generated neutron momentum for Edep ONLY
		double p_n_gen = neutron.Mag(); 
		nEdep = fBAND->GetEdep(p_n_gen);

		// Now use neutron angles and TOF to get neutron variables
		theta_n = neutron.Theta()*TMath::RadToDeg();
		phi_n = neutron.Phi()*TMath::RadToDeg();
		double dZ_n = -fRand->Uniform(zDown + vz, zUp + vz);  
		dL_n = dZ_n/cos(neutron.Theta()); 	
		nTime = fRand->Gaus(fRand->Uniform(-60., 330.), 8.5);
		nBeta = (dL_n/100.)/(c*nTime*1.e-9);
		p_n = Mn / sqrt( 1./pow(nBeta,2) - 1. );
		E_n = sqrt(Mn*Mn + p_n*p_n );

		// Photon variables
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
		phi_nq   = norm_scatter.Angle( norm_reaction );
		theta_nq = neutron.Angle( q );
		CosTheta_nq = cos(theta_nq);
		TVector3 direction = norm_scatter.Cross(norm_reaction);
	
		if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
		}
		else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
			phi_nq *= (-1);
		}

		// Remaining kinematic variables
		xB = Q2/(2*Mp*nu);
		W2 =  Mp*Mp + 2.*nu*Mp - Q2;
		double W_primeSq = MD*MD - Q2 + Mn*Mn + 2.*MD*(nu-E_n) - 2.*nu*E_n + 2.*q.Mag()*p_n*cos(theta_nq);
		Wp = sqrt(W_primeSq);
		Xp = Q2/(2.*( nu*(MD-E_n) + p_n*q.Mag()*CosTheta_nq));
		As = (E_n - p_n*CosTheta_nq)/Mn;

		// Check acceptance for electron, neutron
		sector_e = fCLAS->GetElectronAcceptance(theta_e, phi_e, p_e);
		double acc_n = fBAND->GetNeutronAcceptance(theta_n, phi_n, p_n, vz);
		double p2b = fBAND->PointsToBAND(theta_n, phi_n, vz);

		if (sector_e > 0 && p2b > 0.5) {
			
			SetEventData();
			fIO->FillTree();
			nBG++;

		}

	return nBG;

}

void Simulate::SetEventData() {

	fIO->fCLASHit.setSector(sector_e);
	fIO->fCLASHit.setPID(11);
	fIO->fCLASHit.setCharge(-1.);
	fIO->fCLASHit.setStatus(1);
	fIO->fCLASHit.setTime(0.);
	fIO->fCLASHit.setBeta(eBeta);
	fIO->fCLASHit.setChi2(1.);
	fIO->fCLASHit.setEtot(E_e);
	fIO->fCLASHit.setEpcal(E_e);
	fIO->fCLASHit.setEoP(1);
	fIO->fCLASHit.setTimeScint(1);
	fIO->fCLASHit.setPathScint(1);
	fIO->fCLASHit.setU(1);
	fIO->fCLASHit.setV(1);
	fIO->fCLASHit.setQ(1);
	fIO->fCLASHit.setVtx(vx);
	fIO->fCLASHit.setVty(vy);
	fIO->fCLASHit.setVtz(vz);
	fIO->fCLASHit.setMomentum(p_e);
	fIO->fCLASHit.setTheta(theta_e);
	fIO->fCLASHit.setPhi(phi_e);
	fIO->fCLASHit.setQ(0);
	fIO->fCLASHit.setThetaQ(0);
	fIO->fCLASHit.setPhiQ(0);
	fIO->fCLASHit.setQ2(Q2);
	fIO->fCLASHit.setOmega(nu);
	fIO->fCLASHit.setXb(xB);
	fIO->fCLASHit.setW2(W2);

	fIO->fBANDHit.setSector(0);
	fIO->fBANDHit.setLayer(0);
	fIO->fBANDHit.setComponent(0);
	fIO->fBANDHit.setBarID(0);
	fIO->fBANDHit.setEdep(nEdep);
	fIO->fBANDHit.setTof(nTime);
	fIO->fBANDHit.setTofFadc(nTime);
	fIO->fBANDHit.setTdiff(nTime);
	fIO->fBANDHit.setTdiffFadc(nTime);
	fIO->fBANDHit.setX(0);
	fIO->fBANDHit.setY(0);
	fIO->fBANDHit.setZ(0);
	fIO->fBANDHit.setStatus(1);
	fIO->fBANDHit.setRawLtdc(1);		
	fIO->fBANDHit.setRawRtdc(1);		
	fIO->fBANDHit.setRawLtdccorr(1);	
	fIO->fBANDHit.setRawRtdccorr(1);	
	fIO->fBANDHit.setRawLtfadc(1);	
	fIO->fBANDHit.setRawRtfadc(1);	
	fIO->fBANDHit.setRawLamp(1);		
	fIO->fBANDHit.setRawRamp(1);		
	fIO->fBANDHit.setRawLadc(1);		
	fIO->fBANDHit.setRawRadc(1);		
	fIO->fBANDHit.setPmtLtdc(1);		
	fIO->fBANDHit.setPmtRtdc(1);		
	fIO->fBANDHit.setPmtLtfadc(1);	
	fIO->fBANDHit.setPmtRtfadc(1);	
	fIO->fBANDHit.setPmtLamp(1);		
	fIO->fBANDHit.setPmtRamp(1);		
	fIO->fBANDHit.setPmtLadc(1);		
	fIO->fBANDHit.setPmtRadc(1);		
	fIO->fBANDHit.setPmtLped(1);		
	fIO->fBANDHit.setPmtRped(1);		

	fIO->fTagHit.setThetaNQ(theta_nq);
	fIO->fTagHit.setPhiNQ(phi_nq);
	fIO->fTagHit.setWp(Wp);
	fIO->fTagHit.setXp(Xp);
	fIO->fTagHit.setAs(As);

	fIO->bg = bg;

}

