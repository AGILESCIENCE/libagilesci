


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




double Sinaa(double angle)
{
return angle ? sin(angle)/angle : 0.0;
}


double SphDistDeg(double long1, double lat1, double long2, double lat2)
{
    double l1 = long1*DEG2RAD;
    double l2 = long2*DEG2RAD;
    double b1 = lat1*DEG2RAD;
    double b2 = lat2*DEG2RAD;
    double dist = acos(sin(b1)*sin(b2) + cos(b1)*cos(b2) *cos(l1-l2));
    return dist*RAD2DEG;
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
