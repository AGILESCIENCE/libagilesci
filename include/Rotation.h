////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       2009-06-03
//       Authors: Tomaso Contessi, IASF-Milano
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

#ifndef _ROTATION_
#define _ROTATION_

#include <cmath>


const double DegToRad = M_PI/180;
const double RadToDeg = 180/M_PI;


class Vect;
class Rotation;


/// /////////////////////////////////////////////
///
/// Class Vect
///


class Vect
{
public:
	double x;
	double y;
	double z;
public:
	Vect() { x = y = z = 0; }
	Vect(double i, double j, double k) { x=i; y=j; z=k; }
	Vect(const double* array) { x=array[0]; y=array[1]; z=array[2]; }
	Vect(const Vect& v) { x=v.x; y=v.y; z=v.z; }

	Vect& operator=(const Vect& v) { x=v.x; y=v.y; z=v.z; return *this; }

	double Norm() const { return sqrt(x*x+y*y+z*z); }

	/// Scaling a vector
	friend Vect operator*(double a, const Vect& v) { Vect v2(v); v2.x*=a; v2.y*=a; v2.z*=a; return v2; }
	friend Vect operator*(const Vect& v, double a) { return a*v; }
	Vect& operator*=(double a) { x*=a; y*=a; z*=a; return *this; }

	/// Rotating a vector
	Vect& operator*=(const Rotation& rot);
	Vect& operator/=(const Rotation& rot);

	friend Rotation operator*(const Rotation& r1, const Rotation& r2);
	friend Rotation operator/(const Rotation& r1, const Rotation& r2);
};


/// /////////////////////////////////////////////
///
/// Class Rotation
///
/// Convention for euler angles:
/// phi   = heading in degrees, applied first
/// theta = attitude in degrees, applied sedond
/// psi   = bank in degrees, applied third
///
/// Conventions:
/// When a vector v of class Vect is used to represent a set of euler angles
/// it is intended in radians, and v.x = heading, v.y = attitude, v.z = bank
///


class Rotation
{
public:
	Rotation() { w=1; x=y=z=0; }
	Rotation(double quat_w, double quat_x, double quat_y, double quat_z) { w = quat_w; x = quat_x; y = quat_y; z = quat_z; }
	Rotation(double phi, double theta, double psi) { SetEuler(phi, theta, psi); }		/// Degrees
	Rotation(Vect eulerRad) { SetEuler(eulerRad); }												/// Radians
	Rotation(const Rotation& rot) { w = rot.w; x = rot.x; y = rot.y; z = rot.z; }
	
	Rotation& operator=(const Rotation& rot) { w = rot.w; x = rot.x; y = rot.y; z = rot.z; return *this; }
	
	/// Conversion between Euler angles in Degrees
	void SetEuler(double phi, double theta, double psi);
	void GetEuler(double& phi, double& theta, double& psi) const;

	/// Conversion between Euler angles in Radians
	void SetEuler(Vect euler);
	Vect GetEuler() const;

	/// Inverse rotations
	Rotation& Invert() { w= -w; return *this; }
	Rotation  Inverse() const { return Rotation(-w,x,y,z); }

	/// Rotation operators
	Rotation& operator-() { w = -w; return *this; }
	Rotation& operator*=(const Rotation& rot) { return *this = *this * rot; }
	Rotation& operator/=(const Rotation& rot) { return *this = *this * rot.Inverse(); }

	
	friend Rotation operator*(const Rotation& r1, const Rotation& r2);
	friend Rotation operator/(const Rotation& r1, const Rotation& r2) { return r1 * r2.Inverse(); }
	friend Vect operator*(const Vect& v, const Rotation& rot);
	friend Vect operator/(const Vect& v, const Rotation& rot) { return v * rot.Inverse(); }

public:
	double w;
	double x;
	double y;
	double z;
};



#endif
