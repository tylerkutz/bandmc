#ifndef IO_HH
#define IO_HH

class TString;
class TFile;
class TTree;

class IO{

public:
	IO(TString);
	~IO();
	
	void FillTree();
	void WriteTree();

	double p_e;
	double theta_e;
	double phi_e;
	double p_n;
	double theta_n;
	double phi_n;


private:

	TFile* fFile;
	TTree* fTree;

	void InitializeTree();


};

#endif

