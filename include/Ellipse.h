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
	Ellipse(const Ellipse& ellipse) { horAxis=ellipse.horAxis; verAxis=ellipse.verAxis; attitude=ellipse.attitude; center=ellipse.center;}
	~Ellipse() {}
	Ellipse& operator=(const Ellipse& ellipse) { horAxis=ellipse.horAxis; verAxis=ellipse.verAxis; attitude=ellipse.attitude; center=ellipse.center; return *this; }
};

/// Find the best fitting ellipse starting from a circle of the same area and centered in the barycente
int PolyToEllipse(const Polygon& p, Ellipse& e);
int PolyToEllipse(const Polygon& p, double area, const Point& barycenter, Ellipse& e);

/// Find the best fitting ellipse starting from a circle of given center and radius
int PolyToEllipse(const Polygon& p, const Point& center, double radius, Ellipse& e);


#endif
