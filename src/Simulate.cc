#include "Simulate.hh"
#include "IO.hh"
#include "CLAS.hh"
#include "BAND.hh"

#include "clashit.h"
#include "bandhit.h"
#include "taghit.h"

#include "TVector3.h"
#include "TRandom3.h"


#include <iostream>
using namespace std;

Simulate::Simulate(bool do_smear, bool inclusive = false) {

	doSmearing = do_smear;
	isInclusive = inclusive;

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

int Simulate::SimulateEvent(TVector3 electron, TVector3 neutron, double radiation) {

		fIO->ClearEvent();
		
		bg = 0;
		radweight = radiation;

		// Vertex
		vx = fRand->Uniform(targetR) * cos(fRand->Uniform(360.*TMath::DegToRad()));
		vy = fRand->Uniform(targetR) * sin(fRand->Uniform(360.*TMath::DegToRad()));
		vz = fRand->Uniform(-8., 3.);

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
		X_n = dL_n * sin(neutron.Theta() ) * cos(neutron.Phi() );
		Y_n = dL_n * sin(neutron.Theta() ) * sin(neutron.Phi() );
		Z_n = dL_n * cos(neutron.Theta() );

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

		// Momentum vectors
		momentumE = e4.Vect();
		momentumN = neutron;
		momentumQ = q;

		// Check acceptance for electron, neutron
		sector_e = fCLAS->GetElectronAcceptance(theta_e, phi_e, p_e);
		double acc_n = fBAND->GetNeutronAcceptance(theta_n, phi_n, p_n, vz);

		if (sector_e > 0 && (acc_n > 0.5 || isInclusive) ) {

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
		vz = fRand->Uniform(-8., 3.);

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
		X_n = dL_n * sin(neutron.Theta() ) * cos(neutron.Phi() );
		Y_n = dL_n * sin(neutron.Theta() ) * sin(neutron.Phi() );
		Z_n = dL_n * cos(neutron.Theta() );

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
	
		// Momentum vectors
		momentumE = e4.Vect();
		momentumN = neutron;
		momentumQ = q;

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

	//////////////
	// Electron //
	//////////////

	// Used in analysis: simulated values
	fIO->fCLASHit.setPID(11);		
	fIO->fCLASHit.setCharge(-1);
	fIO->fCLASHit.setEpcal(E_e);
	fIO->fCLASHit.setVtz(vz);
	fIO->fCLASHit.setQ2(Q2);

	// Used in analysis: dummy values
	fIO->fCLASHit.setEoP(0.2);
	fIO->fCLASHit.setU(16);
	fIO->fCLASHit.setV(16);
	fIO->fCLASHit.setW(16);
	
	// Not used in analysis
	fIO->fCLASHit.setSector(sector_e);
	fIO->fCLASHit.setStatus(1);
	fIO->fCLASHit.setTime(0.);
	fIO->fCLASHit.setBeta(eBeta);
	fIO->fCLASHit.setChi2(1.);
	fIO->fCLASHit.setEtot(E_e);
	fIO->fCLASHit.setTimeScint(1);
	fIO->fCLASHit.setPathScint(1);
	fIO->fCLASHit.setVtx(vx);
	fIO->fCLASHit.setVty(vy);
	fIO->fCLASHit.setMomentum(p_e);
	fIO->fCLASHit.setTheta(theta_e*TMath::DegToRad());
	fIO->fCLASHit.setPhi(phi_e*TMath::DegToRad());
	fIO->fCLASHit.setQ(0);
	fIO->fCLASHit.setThetaQ(0);
	fIO->fCLASHit.setPhiQ(0);
	fIO->fCLASHit.setOmega(nu);
	fIO->fCLASHit.setXb(xB);
	fIO->fCLASHit.setW2(W2);

	/////////////
	// Neutron //
	/////////////

	bandhit thisBANDHit;

	// Used in analysis: simulated values
	thisBANDHit.setTof(nTime);
	thisBANDHit.setX(X_n);
	thisBANDHit.setY(Y_n);
	thisBANDHit.setZ(Z_n);

	// Used in analysis: dummy values
	fIO->fnMult = 1;
	thisBANDHit.setStatus(0);
	thisBANDHit.setEdep(12000.); 
	
	// Not used in analysis
	thisBANDHit.setSector(0);
	thisBANDHit.setLayer(0);
	thisBANDHit.setComponent(0);
	thisBANDHit.setBarID(0);
	thisBANDHit.setTofFadc(nTime);
	thisBANDHit.setTdiff(nTime);
	thisBANDHit.setTdiffFadc(nTime);
	thisBANDHit.setRawLtdc(1);		
	thisBANDHit.setRawRtdc(1);		
	thisBANDHit.setRawLtdccorr(1);	
	thisBANDHit.setRawRtdccorr(1);	
	thisBANDHit.setRawLtfadc(1);	
	thisBANDHit.setRawRtfadc(1);	
	thisBANDHit.setRawLamp(1);		
	thisBANDHit.setRawRamp(1);		
	thisBANDHit.setRawLadc(1);		
	thisBANDHit.setRawRadc(1);		
	thisBANDHit.setPmtLtdc(1);		
	thisBANDHit.setPmtRtdc(1);		
	thisBANDHit.setPmtLtfadc(1);	
	thisBANDHit.setPmtRtfadc(1);	
	thisBANDHit.setPmtLamp(1);		
	thisBANDHit.setPmtRamp(1);		
	thisBANDHit.setPmtLadc(1);		
	thisBANDHit.setPmtRadc(1);		
	thisBANDHit.setPmtLped(1);		
	thisBANDHit.setPmtRped(1);		

	new(fIO->saveHit[0]) bandhit;
	fIO->saveHit[0] = &thisBANDHit;

	/////////
	// Tag //
	/////////
	
	taghit thisTagHit;

	thisTagHit.setThetaNQ(theta_nq);
	thisTagHit.setPhiNQ(phi_nq);
	thisTagHit.setWp(Wp);
	thisTagHit.setXp(Xp);
	thisTagHit.setAs(As);

	thisTagHit.setMomentumE(momentumE);
	thisTagHit.setMomentumN(momentumN);
	thisTagHit.setMomentumQ(momentumQ);

	new(fIO->saveTag[0]) taghit;
	fIO->saveTag[0] = &thisTagHit;

	// Misc.
	fIO->bg = bg;
	fIO->radweight = radweight;

}

