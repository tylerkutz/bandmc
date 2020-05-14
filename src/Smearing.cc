#include "Smearing.hh"

#include "TRandom3.h"
#include "TLorentzVector.h"

#include <unistd.h>

using namespace std;

Smearing::Smearing() {

	fRand = new TRandom3(0);
	int seed = time(0) - getpid();
	fRand->SetSeed(seed);

}

Smearing::~Smearing() {
}

void Smearing::smear(TLorentzVector *V4, int q){
        double inM = V4->M();
        double sP  = V4->P();
        double sTh = V4->Theta();
        double sPh = V4->Phi();

        double sThD = TMath::RadToDeg()*sTh;
        double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
        double momS2 = 0.02*sThD;
        double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*sP,2) + TMath::Power(momS2,2) );
        momR *= 2.0;

        double theS1 = 0.004*sThD + 0.1;
        double theS2 = 0;
        double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(theS2,2) );
        theR *= 2.5;

        double phiS1 = 0.85-0.015*sThD;
        double phiS2 = 0.17-0.003*sThD;
        double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(phiS2,2) );
        phiR *= 3.5;

        sPh += TMath::DegToRad() * phiR * fRand->Gaus();
        sTh += TMath::DegToRad() * theR * fRand->Gaus();
        sP  += momR  * fRand->Gaus() *  V4->P() ;

        V4->SetE( TMath::Sqrt( sP*sP + inM*inM )  );
        V4->SetRho( sP );
        V4->SetTheta( sTh );
        V4->SetPhi( sPh );
}

void Smearing::smearCD(TLorentzVector *V4, int q){
        double inM = V4->M();
        double sP  = V4->P();
        double sTh = V4->Theta();
        double sPh = V4->Phi();

        double sThD = TMath::RadToDeg()*sTh;
        double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
        double momS2 = 0.02*sThD;
        double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*sP,2) + TMath::Power(momS2,2) );
        momR *= 5;

        double theS1 = 0.004*sThD + 0.1;
        double theS2 = 0;
        double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(theS2,2) );
        theR *= 1.0;

        double phiS1 = 0.85-0.015*sThD;
        double phiS2 = 0.17-0.003*sThD;
        double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(phiS2,2) );
        phiR *= 3.5;

        sPh += TMath::DegToRad() * phiR * fRand->Gaus();
        sTh += TMath::DegToRad() * theR * fRand->Gaus();
        sP  += momR  * fRand->Gaus() *  V4->P() ;

        V4->SetE( TMath::Sqrt( sP*sP + inM*inM )  );
        V4->SetRho( sP );
        V4->SetTheta( sTh );
        V4->SetPhi( sPh );
}
