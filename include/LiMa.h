////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       Authors: Andrea Bulgarelli, Leonardo Baroncelli, Giancarlo Zollino (INAF/OAS Bologna)
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       Copyright (C) 2005-2019 AGILE Team. All rights reserved.
/*
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
////////////////////////////////////////////////////////////////////////////////////
 
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
