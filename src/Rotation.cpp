////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       2009
//       Authors: Tomaso Contessi (NUOVE Idee SAS)
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



#include "Rotation.h"


/// /////////////////////////////////////////////
///
/// Class Vect Definition
///

Vect& Vect::operator*=(const Rotation& rot)
{
*this = *this * rot;
return *this;
}

Vect& Vect::operator/=(const Rotation& rot)
{
*this = *this * rot.Inverse();
return *this;
}

Vect operator*(const Vect& v, const Rotation& rot)
{
return Vect(
v.x - v.x*(rot.y*rot.y+rot.z*rot.z)*2 + v.y*(rot.x*rot.y-rot.z*rot.w)*2 + v.z*(rot.x*rot.z+rot.y*rot.w)*2,
v.x*(rot.x*rot.y+rot.z*rot.w)*2 + v.y - v.y*(rot.x*rot.x+rot.z*rot.z)*2 + v.z*(rot.y*rot.z-rot.x*rot.w)*2,
v.x*(rot.x*rot.z-rot.y*rot.w)*2 + v.y*(rot.y*rot.z+rot.x*rot.w)*2 + v.z - v.z*(rot.x*rot.x+rot.y*rot.y)*2
);
}


/// /////////////////////////////////////////////
///
/// Class Rotation Definition
///


void Rotation::SetEuler(double phi, double theta, double psi)
{
Vect euler(phi*DegToRad, theta*DegToRad, psi*DegToRad);
SetEuler(euler);
}


/// See http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/index.htm
void Rotation::SetEuler(Vect euler)
{
double cosf = cos(euler.x/2);
double cost = cos(euler.y/2);
double cosp = cos(euler.z/2);
double sinf = sin(euler.x/2);
double sint = sin(euler.y/2);
double sinp = sin(euler.z/2);

w = cosf*cost*cosp - sinf*sint*sinp;
x = sinf*sint*cosp + cosf*cost*sinp;
y = sinf*cost*cosp + cosf*sint*sinp;
z = cosf*sint*cosp - sinf*cost*sinp;
}



void Rotation::GetEuler(double& phi, double& theta, double& psi) const
{
Vect v = GetEuler();
phi = v.x*RadToDeg;
theta = v.y*RadToDeg;
psi = v.z*RadToDeg;
}


/// See http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToEuler/
Vect Rotation::GetEuler() const
{
Vect eul;
double test = 2*(x*y + z*w);
if (test==1.0) {
	eul.x = 2 * atan2(x,w);
	eul.y = M_PI/2;
	eul.z = 0;
	}
else if (test==-1.0) {
	eul.x = -2 * atan2(x,w);
	eul.y = -M_PI/2;
	eul.z = 0;
	}
else {
	eul.x = atan2(2*y*w-2*x*z , 1 - 2*y*y - 2*z*z);
	eul.y = asin(2*x*y + 2*z*w);
	eul.z = atan2(2*x*w-2*y*z , 1 - 2*x*x - 2*z*z);
	}
return eul;
}


Rotation operator*(const Rotation& r1, const Rotation& r2)
{
Rotation q;
q.x =  r1.w*r2.x + r1.z*r2.y - r1.y*r2.z + r1.x*r2.w;
q.y = -r1.z*r2.x + r1.w*r2.y + r1.x*r2.z + r1.y*r2.w;
q.z =  r1.y*r2.x - r1.x*r2.y + r1.w*r2.z + r1.z*r2.w;
q.w = -r1.x*r2.x - r1.y*r2.y - r1.z*r2.z + r1.w*r2.w;
return q;
}
