


#ifndef _MATH_UTILS_
#define _MATH_UTILS_


#include <cmath>

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


double Sinaa(double angle);

double SphDistDeg(double long1, double lat1, double long2, double lat2);

void Euler(double ai, double bi, double * ao, double * bo, int select);

#endif
