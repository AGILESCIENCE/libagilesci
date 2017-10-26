/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/
#include "LiMa.h"


LiMa::LiMa(double _ctsBinSumT0,double _ctsBinSumT1,double _ctsBinSumT2,double _expBinSumT0,double _expBinSumT1,double _expBinSumT2)
{

	ctsBinSumT0 = _ctsBinSumT0;
	ctsBinSumT1 = _ctsBinSumT1;
	ctsBinSumT2 = _ctsBinSumT2;
	expBinSumT0 = _expBinSumT0;
	expBinSumT1 = _expBinSumT1;
	expBinSumT2 = _expBinSumT2;
	

	bkg = ctsBinSumT1 + ctsBinSumT2;
	source = ctsBinSumT0;
	//N_on = source + bkg;
	N_off = bkg;
	expBgSum = expBinSumT1 + expBinSumT2;
	alpha = expBinSumT0 / (expBgSum);
	
	alp1 = alpha / (1 + alpha);
	alp2 = alpha + 1;

	cout << "sig cts (N_on)  " << source << endl;
	cout << "sig exp  " << expBinSumT0 << endl;
	cout << "sig rate " << setprecision(10) << source / (double) expBinSumT0 << endl;
	cout << "bkg cts (N_off) " << bkg << endl;
	cout << "bkg exp  " << (expBgSum) << endl;
	cout << "bkg rate " << setprecision(10) << bkg / (double) (expBgSum) << endl;
	cout << "alpha: " << std::setprecision(10) << alpha << endl;
	cout << "N_s = N_on - alpha * N_off: " << std::setprecision(10) << source - alpha * bkg << endl;
}


double LiMa::computeLiMiValue()
{
 
 
	if ((source > 0) and (bkg > 0)) {
		double source1 = source;
		double bkg1 = bkg;
		double L1 = pow(((source1 + bkg1) / source1) * alp1, source);
		double L2 = pow(((bkg1 + source1) / bkg1) / alp2, bkg);
		double L = L1 * L2;
		S = sqrt(-2. * log(L));
		SA = sqrt(2.) * sqrt(source1 * log( (1 / alp1 ) * ( source1 / (double)(source1 + bkg1) )) + bkg1 * log( alp2 * ( bkg1 / (double)( source1 + bkg1 ) ) ) );
		cout <<  "Li&Ma sigma " << S << endl;
		cout <<  "Li&Ma sigma " << SA << endl;
		
	} else {
		cout << "Alpha: 0" << endl;
		cout << "Li&Ma sigma 0" << endl;
	}
	return S;
}
