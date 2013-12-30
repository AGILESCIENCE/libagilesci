


#ifndef _ELLIPSE_
#define _ELLIPSE_


#include "Polygon.h"


struct Ellipse
{
	double horAxis;
	double verAxis;
	double attitude;	///  radians;
	Point  center;
	Ellipse(): horAxis(0), verAxis(0), attitude(0), center() {}
};

/// Find the best fitting ellipse starting from a circle of the same area and centered in the barycente
int PolyToEllipse(const Polygon& p, Ellipse& e);
int PolyToEllipse(const Polygon& p, double area, const Point& barycenter, Ellipse& e);

/// Find the best fitting ellipse starting from a circle of given center and radius
int PolyToEllipse(const Polygon& p, const Point& center, double radius, Ellipse& e);


#endif
