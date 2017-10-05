/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/


 
 
#include <string>
#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

class LiMa
{
public:
 	LiMa(double ctsBinSumT0,double ctsBinSumT1,double ctsBinSumT2,double expBinSumT0,double expBinSumT1,double expBinSumT2);
	double computeLiMiValue();
 

	double ctsBinSumT0;
	double ctsBinSumT1;
	double ctsBinSumT2;
	double expBinSumT0;
	double expBinSumT1;
	double expBinSumT2;
	
	int bkg;
	int source;
	int N_on;
	int N_off;	
	double expBgSum;
	double alpha;
	double alp1;
	double alp2;
	
	double S;
	double SA;

};
