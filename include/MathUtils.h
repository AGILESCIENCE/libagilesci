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

#ifndef _MATH_UTILS_
#define _MATH_UTILS_


#include <cmath>
#include <assert.h>
#include "TMath.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


const double PI_BY_TWO = M_PI/2;
const double TWO_PI = M_PI*2;
const double FOUR_PI = M_PI*4;

const double RAD2DEG = 180.0/M_PI;
const double DEG2RAD = M_PI/180.0;


/// Convert a double to a nice string. Strip the trailing zeroes by default.
/// Return (and set) the user string or an internal static
char* DoubleStr(double dVal, int decimals=-1, char* dStr=0);

inline double Sinaa(double angle)
{
    return angle ? sin(angle)/angle : 0.0;
}

/// distance in Galactic coordinates. Longitude=l, latitude=b
inline double SphDistDeg(double long1, double lat1, double long2, double lat2)
{
    double l1 = long1*DEG2RAD;
    double l2 = long2*DEG2RAD;
    double b1 = lat1*DEG2RAD;
    double b2 = lat2*DEG2RAD;
    double val = sin(b1)*sin(b2) + cos(b1)*cos(b2) *cos(l1-l2);
    if(val > 1.0) return 0.0;
    if(val < -1.0) return 180.0;
    double dist = acos(val);
    return dist*RAD2DEG;
}

/// convert from Equatorial (RA, DEC) to Galactic (l,b) coordinates
void Euler(double ai, double bi, double * ao, double * bo, int select);

inline double Exposure_cm2s(double exp, double lat, double mres)
{
	//exp [cm2 s sr]
	//return in [cm2 s]
	return exp  / ( pow(mres/180.0,2) * cos(lat));	
}

/// The return value is 0 on success, !=0 on failure.
/// n = number of data points
/// x,y  = arrays of data
/// b = output intercept
/// m  = output slope
/// r = output correlation coefficient (can be NULL if you don't want it)
/// y = mx + b
int linreg(int n, const double x[], const double y[], double& b, double& m, double& r);

inline double CalcLogBarycenter(double e_min, double e_max) {
	double e_log = (TMath::Log10(e_max) - TMath::Log10(e_min)) / 2.0 + TMath::Log10(e_min);
	return pow(10, e_log);
}
#endif
