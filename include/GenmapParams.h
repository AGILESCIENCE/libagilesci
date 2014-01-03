


#ifndef _GENMAP_PARAMS_
#define _GENMAP_PARAMS_


#include <vector>
#include <string>
#include "fitsio.h"
#include "MathUtils.h"
#include "PilParams.h"
#include "Intervals.h"



const double obtlimit = 104407200.0;
enum ProjType {ARC, AIT};


struct MaplistEntry
{
	double fovradmin;
	double fovradmax;
	double emin;
	double emax;
	double index;
};

class Maplist : public std::vector<MaplistEntry>
{
public:
	Maplist() {}
	bool Store(const char* maplist);
	void Print() ;//const;
};



/// EXP


class ExpGenCommonParams: public PilParams
{
public:
	ExpGenCommonParams(const PilDescription params[]): PilParams(params) {}

	bool AssignCommonParams();

	bool YTolTest(double ra_y, double dec_y, double ra_y0, double dec_y0) const
		{ 	return SphDistDeg(ra_y, dec_y, ra_y0, dec_y0) > y_tol; }
	bool EarthTolTest(double earth_ra, double earth_dec, double earth_ra0, double earth_dec0) const
		{ return SphDistDeg(earth_ra, earth_dec, earth_ra0, earth_dec0) > earth_tol; }
	bool FovTest(double theta) const { return (theta < fovradmax && theta >= fovradmin); }
	bool AlbTest(double ra, double dec, double earth_ra, double earth_dec) const
		{ return SphDistDeg(ra, dec, earth_ra, earth_dec) > albrad; }
	bool RollTolTest(double * psi) {return fabs(psi[0] - psi[1]) > roll_tol;}

	Intervals intervals;
	const char* outfile;
	const char* logfile;
	const char* sarFileName;
	const char* timelist;
	double mdim;
	double la;
	double ba;
	double lonpole;
	double albrad;
	double y_tol;
	double roll_tol;
	double earth_tol;
	int phasecode;
	ProjType projection;
	double timestep;
	double index;
	double tmin;
	double tmax;
	double emin;
	double emax;
	double fovradmin;
	double fovradmax;
};


class ExpGenParams: public ExpGenCommonParams
{
public:
	ExpGenParams();
	bool Load(int argC, char* argV[]);

	int write_fits_header(long i, fitsfile* mapFits, int& status);
	/// bool inmap(int i, int ii) {return (i < mxdim) && (i >= 0) && (ii < mxdim) && (ii >= 0);}
	/// bool etest(long i, double e) {return (e >= maps[i].emin) && (e <= maps[i].emax);}
	bool FovTest(long i, double theta) { return theta < maps[i].fovradmax && theta >= maps[i].fovradmin; }
public:
	Maplist   maps;
	const char* edpFileName;
	const char* maplist;
	double mres;
	double binstep;
	long mxdim;
	};


class ExpGenParamsT: public ExpGenCommonParams
{
public:
	ExpGenParamsT();
	bool Load(int argC, char* argV[]);
public:
	double timeslot;
	};

class GammaExtractParams: public ExpGenCommonParams
{
public:
	GammaExtractParams();
	bool Load(int argC, char* argV[]);
public:
	double timeslot;
	const char* evtfile;
	int filtercode;
	};
	
/// CTS


class CtsGenCommonParams: public PilParams
{
public:
	CtsGenCommonParams(const PilDescription params[]): PilParams(params) {}
	bool AssignCommonParams();
	bool inmap(int i, int ii) { return (i < mxdim) && (i >= 0) && (ii < mxdim) && (ii >= 0); }
	int write_fits_header(fitsfile *mapFits, ProjType projection, int & status);
public:
	Intervals intervals;
	const char* evtfile;
	const char* outfile;
	const char* timelist;
	double mdim;
	double mres;
	long mxdim;
	double la;
	double ba;
	double tmin;
	double tmax;
	double emin;
	double emax;
	double fovradmax;
	double albrad;
	double lonpole;
	int phasecode;
	int filtercode;
};

class CtsGenParams: public CtsGenCommonParams
{
public:
	CtsGenParams();
	bool Load(int argC, char* argV[]);
public:
	double fovradmin;
	ProjType projection;
};

class CtsGenParamsT: public CtsGenCommonParams
{
public:
	CtsGenParamsT();
	bool Load(int argC, char* argV[]);
public:
	double fovradmin;
	double timeslot;
};

class ThetaGenParams: public CtsGenCommonParams
{
public:
	ThetaGenParams();
	bool Load(int argC, char* argV[]);
public:
	ProjType projection;
};


#endif
