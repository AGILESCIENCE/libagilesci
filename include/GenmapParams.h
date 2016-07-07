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

class ThetaGenParams: public CtsGenCommonParams
{
public:
	ThetaGenParams();
	bool Load(int argC, char* argV[]);
public:
	ProjType projection;
};

#endif
