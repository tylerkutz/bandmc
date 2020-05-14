#include "RadGen.hh"

using namespace std;

extern "C" {

	void radgen_init_(string* CTARGET, float* ebeam, int* LST40, int* ixytst);
	void RADGEN_(float* e1, double* vpgen, double* vprad, float* phrad, float* q2tr, float*anutr, float* weight);

}



RadGen::RadGen() {

	eBeam = 10.6; // GeV
	target = "D";

	radgen_init_(&target, &eBeam, 0, -1);


}

RadGen::~RadGen() {

}

void RadGen::Radiate(double* qvec) {



}


