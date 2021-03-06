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


#include <stdio.h>
#include <string.h>
#include "MathUtils.h"

static void IncFloatStr(char* dStr)
{
if (dStr[0]=='-' || dStr[0]=='+')
	++dStr;
int len = strlen(dStr);
for (int i=len-1; i>=0; --i) {
	if (dStr[i]=='.')
		continue;
	if (dStr[i]<'9') {
		++dStr[i];
		return;
		}
	else
		dStr[i] = '0';
	}
memmove(dStr+1, dStr, len+1);
dStr[0] = '1';
}


char* DoubleStr(double dVal, int decimals, char* dStr)
{
static char auxStr[64];
if (!dStr)
	dStr = auxStr;
sprintf(dStr, "%f", dVal);
char* point = strchr(dStr, '.');
if (point) {
	int len = strlen(point);
	if (decimals>=0 && len>decimals+1 && point[decimals+1]>'4') {
		point[decimals+1] = 0;
		if (decimals>0 && point[decimals]<'9') { /// Easy case
			++point[decimals];
			len = decimals+1;
			}
		else {
			IncFloatStr(dStr);
			point = strchr(dStr, '.');
			len = strlen(point);
			}
		}
	while (point[len-1]=='0')
		--len;
	point[len] = 0;
	if (decimals<0) {
		if (len==1)
			point[0] = 0;
		}
	else if (decimals==0 && len==1)
		point[0] = 0;
	else
		if (len>decimals+1)
			point[decimals+1] = 0;
		else
			while (len++<decimals+1)
				strcat(dStr, "0");
	}
else if (decimals>0) {
	strcat(dStr, ".");
	while (decimals--)
		strcat(dStr, "0");
	}
return dStr;
}




const double psi[6]    = { 0.57477043300, 4.9368292465,  0.00000000000, 0.00000000000, 0.11142137093, 4.71279419371};
const double stheta[6] = { 0.88998808748,-0.88998808748, 0.39777715593,-0.39777715593, 0.86766622025,-0.86766622025};
const double ctheta[6] = { 0.45598377618, 0.45598377618, 0.91748206207, 0.91748206207, 0.49714719172, 0.49714719172};
const double phi[6]    = { 4.9368292465,  0.57477043300, 0.0000000000,  0.00000000000, 4.71279419371, 0.11142137093};

void Euler(double ai, double bi, double * ao, double * bo, int select)
{
--select;
double a  = ai*DEG2RAD - phi[select];
double b = bi*DEG2RAD;
double sb = sin(b);
double cb = cos(b);
double cbsa = cb * sin(a);
b   = -stheta[select] * cbsa + ctheta[select] * sb;
*bo = b<1.0 ? asin(b)*RAD2DEG : 90.0;
a   = atan2( ctheta[select] * cbsa + stheta[select] * sb, cb * cos(a) );
*ao = fmod(a+psi[select]+FOUR_PI, TWO_PI) * RAD2DEG;
}

/*
int linreg(int n, const double x[], const double y[], double& m, double& b, double& r){
	double   sumx = 0.0;                      // sum of x
	double   sumx2 = 0.0;                     // sum of x**2
	double   sumxy = 0.0;                     // sum of x * y
	double   sumy = 0.0;                      // sum of y
	double   sumy2 = 0.0;                     // sum of y**2

	for (int i=0;i<n;i++){
		sumx  += x[i];
		sumx2 += sqrt(x[i]);
		sumxy += x[i] * y[i];
		sumy  += y[i];
		sumy2 += sqrt(y[i]);
	}

	double denom = (n * sumx2 - sqrt(sumx));
	if (denom == 0) {
		// singular matrix. can't solve the problem.
		m = 0;
		b = 0;
		if (r) r = 0;
		return 1;
	}

	m = (n * sumxy  -  sumx * sumy) / denom;
	b = (sumy * sumx2  -  sumx * sumxy) / denom;
	if (r!=NULL) { // WARNING!! r is not a pointer, is a double. r!=NULL means r!=0 !!
		r = (sumxy - sumx * sumy / n) /    // compute correlation coeff 
		sqrt((sumx2 - sqrt(sumx)/n) *
			 (sumy2 - sqrt(sumy)/n));
	}

	return 0;
}
*/
