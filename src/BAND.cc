#include "BAND.hh"

#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"

#include <math.h>

#include <iostream>

using namespace std;

BAND::BAND() {

	fRand = new TRandom3(0);

	nEffFunc = new TF1("nEffFunc", "2.5*pol6", 0.1, 0.65);
	double funcPars[] = {	-1.33158116e+00,
				2.59593847e+01, 
				-1.81677059e+02,  
				6.77091502e+02, 
				-1.35907420e+03,  
				1.37637068e+03, 
				-5.50389370e+02}; // 5 MeV threshold

	nEffFunc->SetParameters(funcPars);


}

BAND::~BAND(){
}


double BAND::GetNeutronAcceptance(double theta, double phi, double pNeutron, double z_m) {

	double p2B = PointsToBAND(theta, phi, z_m);
	double nEff = Efficiency(pNeutron);

	return p2B*nEff;

}

double BAND::Efficiency(double p) {

	double eff = 0.;

	if(p < 0.1) {
		eff = 0.;
	} else if (p > 0.65) {
		eff = 0.3;
	} else {
		eff = nEffFunc->Eval(p);
	}

	if(fRand->Uniform() < eff) {
		return 1.0;
	} else {
		return 0.0;
	}

}

double BAND::PointsToBAND(double theta, double phi, double z_m) {

        //double z = z_m*100; // from m to cm
	double z = z_m;
	theta*=TMath::DegToRad();
	phi*=TMath::DegToRad();

        // Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
        double thickness  = 7.2;                                // thickness of each bar (cm)
        double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6

        // Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java
        double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

        // Distance from ideal target to upstream end of BAND
        // (from BAND survey report, 02/18/2019)
        double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]

        // Distance from ideal target to downstream end of layer 5
        double zDown = (zUpst + 5*thickness) - z_m;

        double rho   = zDown/cos(theta);
        double xDown = rho*sin(theta)*cos(phi);
        double yDown = rho*sin(theta)*sin(phi);

        double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
        double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

        // Sector boundaries
        double topSec1  = globalY + 13*thickness;
        double topSec2  = globalY + 10*thickness;
        double topSec34 = globalY +  3*thickness;
        double topSec5  = globalY -  3*thickness;
        double downSec5 = globalY -  5*thickness;

        if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

        if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
                        (yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
                        (yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
                        (yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
          ) {
                return 1.0;
	}
 
 	return 0.0;

}



