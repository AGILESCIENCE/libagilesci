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

#ifndef _AGILE_MAP_
#define _AGILE_MAP_

#include "ArrayOf.h"
#include "MathUtils.h"




class AgileMap: public MatD
{
public:	/// Construction
	AgileMap(): MatD(), m_xbin(-1), m_ybin(1), m_x0(1), m_y0(1), m_la2(0.0), m_ba2(0.0)
					, m_lonpole(180.0), m_emin(400), m_emax(50000)
					, m_lp(-999), m_bp(-999), m_gp(0.0), m_mapIndex(2.1)
					, m_fovMin(0.0), m_fovMax(0.0), m_albedo(0.0), m_phaseCode(0l)
					, m_step(0.0), m_tstart(0.0), m_tstop(0.0) {
					    m_dateObs[0]=0; m_dateEnd[0]=0;
					    m_skyL[0]=0; m_skyH[0]=0;
					    m_fileName[0]=0;
					}
	AgileMap(const AgileMap &another): MatD() { Copy(another); }
	AgileMap(const char* fileName): MatD() { Read(fileName); }
	~AgileMap() {}

public:	/// Operators
	AgileMap& operator=(const AgileMap& another) { if (this!=&another) Copy(another); return *this; }

public:	/// I/O (the image is transposed on file)
	int Read(const char* fileName);
	int Write(const char* fileName) const;
	int WriteWithAllMetadata(const char* fileName) const;
	void Zero() { MatD::operator=(0.0); }

public:	/// General information
	bool AlignedTo(const AgileMap& another) const;	/// The two maps perfectly overlap
	bool SameParamsAs(const AgileMap& another) const;	/// AlignedTo and same emin and emax too

	const char* GetFileName() const { return m_fileName; } /// What file was opened

	/// Getting members
	double GetMapCenterL() const { return m_la2; }
	double GetMapCenterB() const { return m_ba2; }
	double GetXbin() const { return m_xbin; }
	double GetYbin() const { return m_ybin; }
	double GetX0() const { return m_x0; }
	double GetY0() const { return m_y0; }

	double GetEmin() const { return m_emin; }
	double GetEmax() const { return m_emax; }
	double GetMapIndex() const { return m_mapIndex; }
	double GetLpoint() const { return m_lp; }
	double GetBpoint() const { return m_bp; }

	double GetFovMin() const { return m_fovMin; }
	double GetFovMax() const { return m_fovMax; }
	double GetAlbedo() const { return m_albedo; }
	long   GetPhaseCode() const { return m_phaseCode; }
	double GetLonpole() const { return m_lonpole; }

	double GetStep() const { return m_step; }	/// exp map only

	const char* GetStartDate() const { return m_dateObs; }
	const char* GetEndDate() const { return m_dateEnd; }

    double GetTstart() const { return m_tstart; }
    double GetTstop() const { return m_tstop; }

	const char* GetSkyL() const { return m_skyL; }
	const char* GetSkyH() const { return m_skyH; }

	/// Setting members
	void SetMapCenter(double l, double b) { m_la2=l, m_ba2=b; }
	void SetXYbin(double xbin, double ybin) { m_xbin=xbin; m_ybin=ybin; }
	void SetXY0(double x0, double y0) { m_x0=x0; m_y0=y0; }

	void SetEnergy(double emin, double emax) { m_emin=emin; m_emax=emax; }
	void SetMapIndex(double mapIndex) { m_mapIndex=mapIndex; }

	void SetLBpoint(double lp, double bp) { m_lp=lp; m_bp=bp; }

	void SetFov(double fovMin, double fovMax) { m_fovMin=fovMin; m_fovMax=fovMax; }
	void SetAlbedo(double albedo) { m_albedo=albedo; }
	void SetPhaseCode(long phaseCode) { m_phaseCode=phaseCode; }

	void SetStep(double step) { m_step=step; }	/// exp map only

	void SetStartDate(const char* dateObs);
	void SetEndDate(const char* dateEnd);

	void SetTT(double tstart, double tstop) { m_tstart = tstart; m_tstop = tstop; }

	/// Other information
	// double x(int i, int j) const { return (i+1-m_x0)*m_xbin; }
	// double y(int i, int j) const { return (j+1-m_y0)*m_ybin; }
	double x(int i) const { return (i+1-m_x0)*m_xbin; }
	double y(int j) const { return (j+1-m_y0)*m_ybin; }

	// double theta(int i, int j) const { double xx = x(i,j), yy = y(i,j); return sqrt(xx*xx+yy*yy); }
	double theta(int i, int j) const { double xx = x(i), yy = y(j); return sqrt(xx*xx+yy*yy); }
	double phi(int i, int j) const;

	double Area(int i, int j) const;

	int i(int bin) const { return MatD::DimOf(bin, 0); }
	int j(int bin) const { return MatD::DimOf(bin, 1); }
	int Bin(int i, int j) const { return MatD::AbsIndex(i,j); }
	bool GetRowCol(double lDeg, double bDeg, int* row, int* col) const;
	void GetCoords(int row, int col, double* lDeg, double* bDeg) const;

	double l(int i, int j) const { return m_lMat(i, j); }
	double b(int i, int j) const { return m_bMat(i, j); }
	double l(int bin) const { return m_lMat(i(bin), j(bin)); }
	double b(int bin) const { return m_bMat(i(bin), j(bin)); }

	double SrcDist(int bin, double lng, double lat) const
		{ return SphDistDeg(lng, lat, l(bin), b(bin)); }

	double SrcDist(int i, int j, double lng, double lat) const
		{ return SphDistDeg(lng, lat, l(i,j), b(i,j)); }
	MatD SrcDist(double lng, double lat) const;
	
	/// Operations
	//Return the sum of the bins inside lng,lat,radius, -1 if the radius is outisde
	double SumBin(double lng, double lat, double radius) const;
	bool IsRadiusInside(double lng, double lat, double radius) const;

private:	/// Data
	/// Mandatory for this map to make sense
	/// CTYPE1 = GLON-ARC
	/// CTYPE2 = GLAT-ARC
	double m_xbin;		/// CDELT1
	double m_ybin;		/// CDELT2
	double m_x0;		/// CRPIX1
	double m_y0;		/// CRPIX2
	double m_la2;		/// CRVAL1
	double m_ba2;		/// CRVAL2
	/// Specific data for Agile maps
	double m_lonpole;	/// LONPOLE
	double m_emin;		/// MINENG
	double m_emax;		/// MAXENG
	double m_lp;		/// SC-Z-LII
	double m_bp;		/// SC-Z-BII
	double m_gp;		/// SC-LONPL
	double m_mapIndex;/// INDEX
	double m_fovMin;	/// FOVMIN
	double m_fovMax;	/// FOV
	double m_albedo;	/// ALBEDO
	long   m_phaseCode;/// PHASECOD
	double m_step;		/// STEP (exposure maps only)
	char m_dateObs[32];	/// DATE-OBS
	char m_dateEnd[32];	/// DATE-END
    double m_tstart;    /// TSTART
    double m_tstop;     /// TSTOP
    char m_skyL[1024];  /// SKYL
    char m_skyH[1024];  /// SKYH
	char m_fileName[1024];

	void Copy(const AgileMap& another);
private:	/// Cash
	void Eval_lb();
	MatD m_lMat;
	MatD m_bMat;
};


#endif
