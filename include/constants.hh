#ifndef CONST_H
#define CONST_H

#define PI 3.14159265 
#define GRIDP 201
#define MASSN 0.93956536 //0.93827231//
#define MASSP 0.93827231
#define MASSn (MASSN+MASSP)/2.
#define MASSPI 0.13957
#define HBARC 0.197328  //GeV*fm
#define ALPHA 0.00729735253
#define I complex<double>(0.,1.)
#define INVHBARC 5.06770453255
#define MUP 2.79
#define DEGRTORAD 0.0174532925
#define RADTODEGR 57.2957795
#define MASSD 1.8756 
#define max(a,b) (((a)>(b))?(a):(b))

// Units of GeV
const double mP = 0.93827231;
const double mN = 0.93956536;
const double mD = 1.8756;
const double Ebeam = 10.6;
const double mE = 0.0005109989;

// Units of cm / ns
const double cAir = 29.9792458;

#endif
