#ifndef SMEARING_HH
#define SMEARING_HH

class TLorentzVector;
class TRandom3;

class Smearing{

public:
	Smearing();
	~Smearing();
	void smear(TLorentzVector *V4, int q); 
	void smearCD(TLorentzVector *V4, int q); 

private:
	TRandom3 *fRand;

};

#endif

