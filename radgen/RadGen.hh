#ifndef RADGEN_HH
#define RADGEN_HH

#include <string>

class RadGen{

public:
	RadGen();
	~RadGen();
	
	void Radiate(double*);

private:
	double eBeam;
	std::string target;

};

#endif

