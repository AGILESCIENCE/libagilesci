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

#include "Polygon.h"

#ifndef M_PI
const double M_PI = 3.14159265358979323846264338327950288;
#endif


void Swap(double& d1, double& d2)
{
double temp = d1;
d1 = d2;
d2 = temp;
}



double Point::DistanceTo(const Point& p) const
{
double s1 = p.x - x;
double s2 = p.y - y;
return sqrt(s1*s1+s2*s2);
}




static Segment::CrossMode Within(double number, double bound1, double bound2)
{
if (number==bound1 || number==bound2)
	return Segment::Edge;
else {
	if (bound2<bound1)
		Swap(bound1, bound2);
	if (number>bound1 && number<bound2)
		return Segment::Full;
	return Segment::None;
	}
}



double Segment::AngleTo(const Point& p) const
{
double angle = b.AngleTo(p)-a.AngleTo(b);
if (angle>M_PI)
	angle -= 2*M_PI;
else if (angle<-M_PI)
	angle += 2*M_PI;
return angle;
}


double Min(double a, double b)
{
return a<b ? a : b;
}

double Max(double a, double b)
{
return a>b ? a : b;
}


void Segment::LineCoeffs(double& A, double& B, double& C) const	/// Ax + By + C = 0
{
A = a.y-b.y;
B = b.x-a.x;
C = -a.y*B-a.x*A;
}


static Segment::CrossMode Intersection(double p1, double p2, double p3, double p4, double& p)
{
if (p1>p2)
	Swap(p1, p2);
if (p3>p4)
	Swap(p3, p4);
if (p1>p3) {
	Swap(p1, p3);
	Swap(p2, p4);
	}
if (p2<p3)
	return Segment::None;
p = p2;
if (p2==p3)
	return Segment::Edge;
return Segment::Full;
}


Segment::CrossMode Segment::Cross(const Segment& other, Point& crossPoint) const
{
double a1, b1, c1;
LineCoeffs(a1, b1, c1);
double a2, b2, c2;
other.LineCoeffs(a2, b2, c2);
double det = a1*b2-a2*b1;
if (det==0) {	/// parallel lines
	if (c1==c2) {	/// segments in the same line
		CrossMode horCross = Intersection(a.x, b.x, other.a.x, other.b.x, crossPoint.x);
		CrossMode verCross = Intersection(a.y, b.y, other.a.y, other.b.y, crossPoint.y);
		if (horCross==None || verCross==None)
			return None;
		if (horCross==Full && verCross==Full)
			return Full;
		return Edge;
		}
	else
		return None;
	}
double x = (b1*c2-b2*c1)/det;
double y = (a2*c1-a1*c2)/det;
CrossMode horCross1 = Within(x, a.x, b.x);
CrossMode verCross1 = Within(y, a.y, b.y);
CrossMode horCross2 = Within(x, other.a.x, other.b.x);
CrossMode verCross2 = Within(y, other.a.y, other.b.y);
if (horCross1==None || verCross1==None || horCross2==None || verCross2==None)
	return None;
crossPoint.x = x;
crossPoint.y = y;
if (horCross1==Full && verCross1==Full && horCross2==Full && verCross2==Full)
	return Full;
return Edge;
}



Polygon::Polygon(const Polygon& other): count(other.count)
{
pointList = new Point[count];
for (int i=0; i<count; ++i)
	pointList[i] = other.pointList[i];
}

Polygon::Polygon(const Point& p1, const Point& p2, const Point& p3): count(3)
{
pointList = new Point[3];
pointList[0] = p1;
pointList[1] = p2;
pointList[2] = p3;
}

Polygon::Polygon(const double* points, int numPoints)
: count(numPoints)
{
pointList = new Point[count];
for (int i=0; i<count; ++i) {
	pointList[i].x = points[i*2];
	pointList[i].y = points[i*2+1];
	}
}

Polygon& Polygon::operator=(const Polygon& other)
{
delete[] pointList;
count = other.count;
pointList = new Point[count];
for (int i=0; i<count; ++i)
	pointList[i] = other.pointList[i];
return *this;
}


bool Polygon::Clockwise() const
{
double angle = 0;
for (int i=0; i<count; ++i) {
	Segment seg(NthPoint(i), NthPoint(i+1));
	angle += seg.AngleTo(NthPoint(i+2));
	}
return angle<0;
}


Polygon::Polygon(const Polygon& other, int exclPoint)
{
while (exclPoint<0)
	exclPoint += count; 
exclPoint %= other.count;
count = other.count-1;
pointList = new Point[count];
int offset = 0;
for (int i=0; i<count; ++i) {
	if (i==exclPoint)
		offset = 1;
	pointList[i] = other.pointList[i+offset];
	}
}



static bool InternalPoint(const Polygon& triangle, bool clockwise,  const Point& p)
{
Polygon t2(3);
t2[2] = p;
for (int i=0; i<3; ++i) {
	t2[0] = triangle[i];
	t2[1] = triangle[i+1];
	if (t2.Clockwise()!=clockwise)
		return false;
	}
return true;
}


/// static int cycleCount;

void Polygon::GetBary(double& area, Point& barycenter) const
{

/// cycleCount = 0;

area = 0;
barycenter = Point();
if (count>=3)
	GetBary(0, area, barycenter);
if (area>0)
	barycenter /= area;
}

void Polygon::GetBary(int point, double& area, Point& barycenter) const
{
if (count<3)
	return;
else if (count==3) {
	double s[3];	// Heron's formula
	Segment med[2];
	for (int i=0; i<3; ++i) {
		Segment seg(NthPoint(i), NthPoint(i+1));
		s[i] = seg.Length();
		if (i<2)
			med[i] = Segment(seg.Mean(), NthPoint(i+2));
		}
	double S = (s[0]+s[1]+s[2])/2;
	double newArea = sqrt(S*(S-s[0])*(S-s[1])*(S-s[2]));
	Point newBarycenter;
	med[0].Cross(med[1], newBarycenter);
/**
	++cycleCount;
	cout << cycleCount << ") " <<
		pointList[0].x << ", " << pointList[0].y << " * " << 
		pointList[1].x << ", " << pointList[1].y << " * " << 
		pointList[2].x << ", " << pointList[2].y << " center: " <<
		newBarycenter.x << ", " << newBarycenter.y << " area: " <<
		newArea << endl;
*/
	area += newArea;
	newBarycenter *= newArea;
	barycenter += newBarycenter;
	}
else {
	bool clockwise = Clockwise();
	bool found = false;
	for (int i=point; i<count+point && !found; ++i) {
		Polygon t(NthPoint(i), NthPoint(i+1), NthPoint(i+2));
		if (clockwise==t.Clockwise()) {
			for (int j = i+3; j<i+count && !found; ++j) {
				if (!InternalPoint(t, clockwise, NthPoint(j))) {
					t.GetBary(0, area, barycenter);
					Polygon p2(*this, i+1);
					p2.GetBary(i, area, barycenter);
					found = true;
					}
				}
			}
		}
	}
}
