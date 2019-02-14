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

#include "Ellipse.h"
#include "TH1D.h"
#include "TF1.h"


const double RadToDeg = 180/M_PI;
const double DegToRad = M_PI/180;


static Point s_center;
static double* s_angleArr = 0;


static Double_t FitEllipseFunc(Double_t *x, Double_t *par)
{
double angle = s_angleArr[int(*x)]-par[2];
// double sint = sin(angle);
// double cost = cos(angle);
// double result = sqrt(par[0]*par[0]*cost*cost + par[1]*par[1]*sint*sint);
// double result = 1.0/(cost*cost/(par[0]*par[0]) + sint*sint/(par[1]*par[1]));
double term1 = cos(angle)/par[0];
double term2 = sin(angle)/par[1];
double result = 1.0/(term1*term1 + term2*term2);
return result; 
}


int PolyToEllipse(const Polygon& p, const Point& center, double radius, Ellipse& e)
{
if (s_angleArr)
	return -1;
s_center = center;
int sideCount = p.Sides();
s_angleArr = new double[sideCount];
TH1D histog("Ellipse","Ellipse", sideCount, -0.5, sideCount-0.5);
double maxDistance = 0;
double minDistance = s_center.DistanceTo(p[0]);
for (int i=0; i<sideCount; i++) {
	double distance = s_center.DistanceTo(p[i]);
	histog.SetBinContent(i+1, distance*distance);
	if (distance>maxDistance)
		maxDistance = distance;
	if (distance<minDistance)
		minDistance = distance;
	s_angleArr[i] = s_center.AngleTo(p[i]);
	}
TF1 func("Model", FitEllipseFunc, 0, 0,  3);
func.SetParName(0, "Axis 1");
func.SetParName(1, "Axis 2");
func.SetParName(2, "Angle");
func.SetParLimits(0, minDistance*0.8, maxDistance*1.2);
func.SetParLimits(1, minDistance*0.8, maxDistance*1.2);
func.SetParLimits(2, -M_PI/2, M_PI/2);
func.SetParameter(0, radius);
func.SetParameter(1, radius);
func.SetParameter(2, 0);
Int_t fitResult = histog.Fit(&func, "NEM" ,"", 0, sideCount);
e.horAxis = func.GetParameter(0);
e.verAxis = func.GetParameter(1);
e.attitude = func.GetParameter(2);
e.center = s_center;
delete[] s_angleArr;
s_angleArr = 0;
return fitResult;	/// zzz test fitResult to give a better result value
}

int PolyToEllipse(const Polygon& p, Ellipse& e)
{
double area;
Point barycenter;
p.GetBary(area, barycenter);
double radius = sqrt(area/M_PI);
return PolyToEllipse(p, barycenter, radius, e);
}
