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

#ifndef _POLYGON_
#define _POLYGON_


#include <cmath>


class Point
{
public:
	double x;
	double y;
public:
	Point(): x(0), y(0) {}
	Point(double point_x, double point_y): x(point_x), y(point_y) {}
	Point(const double* point): x(point[0]), y(point[1]) {}
	Point(const Point& point): x(point.x), y(point.y) {}

	~Point() {}
	Point& operator=(const Point& point) { x=point.x; y=point.y; return *this; }

	Point& operator+=(Point p) { x+=p.x; y+=p.y; return *this; }
	Point& operator*=(double a) { x*=a; y*=a; return *this; }
	Point& operator/=(double a) { x/=a; y/=a; return *this; }

	double DistanceTo(const Point& p) const;
	double AngleTo(const Point& p) const { return atan2(p.y-y, p.x-x); }

	friend Point operator+(const Point& p1, const Point& p2) { return Point(p1.x+p2.x, p1.y+p2.y); }
	friend Point operator-(const Point& p1, const Point& p2) { return Point(p1.x-p2.x, p1.y-p2.y); }
	friend bool operator==(const Point& p1, const Point& p2) { return p1.x==p2.x && p1.y==p2.y; }
	friend bool operator!=(const Point& p1, const Point& p2) { return p1.x!=p2.x || p1.y!=p2.y; }
};



class Segment
{
public:
	Point a;
	Point b;
public:
	Segment(): a(), b() {}
	Segment(const Point& p1, const Point& p2): a(p1), b(p2) {}
	Segment(const Segment& other): a(other.a), b(other.b) {}
	Segment& operator=(const Segment& segment) { a=segment.a; b=segment.b; return *this; }
	double Angle() const { return a.AngleTo(b); }
	double AngleTo(const Point& p) const;
	double Length() const { return a.DistanceTo(b); }

	void LineCoeffs(double& A, double& B, double& C) const;	/// Ax + By + C = 0

	enum CrossMode { None, Edge, Full };
	CrossMode Cross(const Segment& other, Point& crossPoint) const;

	Point Mean() const { return Point((a.x+b.x)/2, (a.y+b.y)/2); }
};



class Polygon
{
public:
	int    count;
	Point* pointList;
public:
	Polygon(): count(0), pointList(0) {}
	Polygon(int numPoints): count(numPoints) { pointList = new Point[numPoints]; }
	Polygon(const double* points, int numPoints);
	Polygon(const Polygon& other);
	Polygon(const Point& p1, const Point& p2, const Point& p3);
	~Polygon() { delete[] pointList; }

	Polygon& operator=(const Polygon& other);

	bool Clockwise() const;
	
	int Sides() const { return count; }
	Segment NthSide(int i) const { return Segment(NthPoint(i), NthPoint(i+1)); }

	Point& NthPoint(int i) { while (i<0) i+=count; return pointList[i%count]; }
	const Point& NthPoint(int i) const { while (i<0) i+=count; return pointList[i%count]; }
	Point& operator[](int i) { return NthPoint(i); }
	const Point& operator[](int i) const { return NthPoint(i); }

	void GetBary(double& area, Point& barycenter) const;
	
private:
	Polygon(const Polygon& other, int exclPoint);
	void GetBary(int point, double& area, Point& barycenter) const;
};


#endif
