////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
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

#ifndef _GENMAP_PARAMS_
#define _GENMAP_PARAMS_


#include <vector>
#include <string>
#include "fitsio.h"
#include "MathUtils.h"
#include "PilParams.h"
#include "Intervals.h"

const double obtlimit = 104407200.0;
enum ProjType {ARC, AIT};

/// CTS

class CtsGenCommonParams: public PilParams
{
public:
	CtsGenCommonParams(const PilDescription params[]): PilParams(params) {}
	bool AssignCommonParams();
	bool inmap(int i, int ii) { return (i < mxdim) && (i >= 0) && (ii < mxdim) && (ii >= 0); }
	int write_fits_header(fitsfile *mapFits, ProjType projection, int & status);
public:
	Intervals intervals;
	const char* evtfile;
	const char* outfile;
	const char* timelist;
	double mdim;
	double mres;
	long mxdim;
	double la;
	double ba;
	double tmin;
	double tmax;
	double emin;
	double emax;
	double fovradmax;
	double albrad;
	double lonpole;
	int phasecode;
	int filtercode;
};

class ThetaGenParams: public CtsGenCommonParams
{
public:
	ThetaGenParams();
	bool Load(int argC, char* argV[]);
public:
	ProjType projection;
};

#endif
