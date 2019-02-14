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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

#include "GenmapParams.h"

using namespace std;

const double mjdrefi = 53005.;
const double mjdreff_OLD = 0.00075444444444444444;
const double mjdreff_2009 = 0.00076601852;
const double mjdreff_2012 = 0.000777593;

static int WriteTime(fitsfile* input, double t1, double t2)
{
int status = 0;
int hdnum = 0;
fits_get_num_hdus(input, &hdnum, &status);

int i;
fits_get_hdu_num(input, &i);

char str[100] = "TT";
float timezero = 0.0;
char secstr[] = "s";

fits_update_key(input, TDOUBLE, "TSTART", &t1,  "[OBT]first event time", &status);
fits_update_key(input, TDOUBLE, "TSTOP", &t2,  "[OBT]last event time", &status);
fits_update_key(input, TSTRING, "TIMESYS",  str, NULL, &status);
fits_update_key(input, TSTRING, "TIMEUNIT",  secstr, NULL, &status);

double mjdreff = mjdreff_OLD;
if (t1 >= 268185600.0)
	mjdreff = mjdreff_2012;
else if (t1 >= 157852800.0)
	mjdreff = mjdreff_2009;

int oldstatus = status;
time_t tt;
time(&tt);
tm* tstr;
tstr = gmtime(&tt);

tstr->tm_year = 107;
tstr->tm_mon = 0;
tstr->tm_mday = 1;
tstr->tm_hour = 0;
tstr->tm_min = 0;
tstr->tm_sec = 0;
tstr->tm_isdst = 0;
tt = timegm(tstr);
double tinc = difftime(tt+1,tt);

tt += time_t(((mjdrefi - 54101 + mjdreff) * 86400 + timezero + t1) / tinc);
tstr = gmtime(&tt);
double inttstart;
double sec = tstr->tm_sec + modf(t1, &inttstart);
char datestr[FLEN_KEYWORD];
fits_time2str( tstr->tm_year + 1900, tstr->tm_mon + 1, tstr->tm_mday, tstr->tm_hour, tstr->tm_min, sec, 0, datestr, &status);
fits_update_key(input, TSTRING, "DATE-OBS",  datestr, "start date and time of the observation(TT)", &status);

tt += (time_t)((t2 - inttstart) / tinc);
tstr = gmtime(&tt);
sec = tstr->tm_sec + modf(t2, &inttstart);
char datestr1[FLEN_KEYWORD];	
fits_time2str( tstr->tm_year + 1900, tstr->tm_mon + 1, tstr->tm_mday, tstr->tm_hour, tstr->tm_min, sec, 0, datestr1, &status);
fits_update_key(input, TSTRING, "DATE-END",  datestr1, "end date and time of the observation(TT)", &status);	

if (i == hdnum) {
	fits_movabs_hdu(input, 1, NULL, &status);		
	fits_update_key(input, TDOUBLE, "TSTART", &t1,  "[OBT]first event time", &status);
	fits_update_key(input, TDOUBLE, "TSTOP", &t2,  "[OBT]last event time", &status);
	fits_update_key(input, TSTRING, "TIMESYS",  str, NULL, &status);
	fits_update_key(input, TSTRING, "TIMEUNIT",  secstr, NULL, &status);
	fits_update_key(input, TFLOAT, "TIMEZERO",  &timezero, str, &status);	
	fits_update_key(input, TSTRING, "DATE-OBS",  datestr, "start date and time of the observation(TT)", &status);
	fits_update_key(input, TSTRING, "DATE-END",  datestr1, "end date and time of the observation(TT)", &status);	
	fits_movabs_hdu(input, i, NULL, &status);		
	}

for (int i=1; i<=hdnum; ++i) {
	fits_movabs_hdu(input, i, NULL, &status);	
	oldstatus = status;
	fits_delete_key(input, "MJDREF", &status);
	if (status == 202)
		status = oldstatus;
	double mjdrefiVal = mjdrefi;
	fits_update_key(input, TDOUBLE, "MJDREFI",  &mjdrefiVal, NULL, &status);
	fits_update_key(input, TDOUBLE, "MJDREFF",  &mjdreff, NULL, &status);
	}
return status;
}

static bool LoadTimeList(const char* timelist, Intervals& intervals, double& tmin, double& tmax)
{
if (strcmp(timelist, "None")) {
	intervals = ReadIntervals(timelist);
	int count = intervals.Count();
	if (!count)
		return false;
	tmin = intervals.Min();
	tmax = intervals.Max();
	cout << "Time intervals:" << endl;
	for (int i=0; i<count; ++i)
		cout << intervals[i].Start() << " " << intervals[i].Stop() << endl;
	}
else {
	Interval intv(tmin, tmax);
	intervals.Add(intv);
	}
if (tmin<=obtlimit) {
	cerr << "ERROR: initial time must be greater than " << obtlimit << endl;
	return false;
	}
return true;
}

static bool LoadProjection(const char* projstr, ProjType& projection)
{
if ((strcmp(projstr,"ARC")==0) || (strcmp(projstr,"arc")==0) ||(strcmp(projstr,"GLON-ARC")==0) || (strcmp(projstr,"GLAT-ARC")==0) )
	projection = ARC;
else if ( (strcmp(projstr,"AIT")==0) || (strcmp(projstr,"ait")==0) ||(strcmp(projstr,"AITOFF")==0) || (strcmp(projstr,"aitoff")==0)  ||(strcmp(projstr,"GLON-AIT")==0) || (strcmp(projstr,"GLAT-AIT")==0) )
	projection = AIT;
else {
	cerr << endl << "ERROR: Invalid Projection" << endl;
	return false;
	}
return true;
}

const PilDescription c_paramsTheta[] = {
	{ PilString, "outfile", "Output file name" },
	{ PilString, "evtfile", "Event file index file name" },
	{ PilReal, "mdim", "Size of Map(degrees)" },
	{ PilReal, "mres", "Bin size(degrees)" },
	{ PilReal, "la", "Longitude of map center (galactic)" },
	{ PilReal, "ba", "Latitude of map center (galactic)" },
	{ PilReal, "lonpole", "Rotation of map (degrees)" },
	{ PilReal, "albrad", "Radius of earth albedo (degrees)" },
	{ PilInt, "phasecode", "Orbital phase code" },
	{ PilInt, "filtercode", "Event filter code" },
	{ PilString, "projection", "Projection (ARC or AIT)" },
	{ PilReal, "tmin", "Initial time(sec)" },
	{ PilReal, "tmax", "Final time(sec)" },
	{ PilReal, "emin", "Min energy" },
	{ PilReal, "emax", "Max energy" },
	{ PilReal, "fovradmax", "Radius of field of view (degrees)" },
	{ PilNone, "", "" }
	};


bool CtsGenCommonParams::AssignCommonParams()
{
PilParams& appParams(*this);

evtfile = appParams["evtfile"];
outfile = appParams["outfile"];
/// timelist = appParams["timelist"];
mdim = appParams["mdim"];
mres = appParams["mres"];
mxdim = long(mdim/mres+0.1);
la = appParams["la"];
ba = appParams["ba"];
tmin = appParams["tmin"];
tmax = appParams["tmax"];
emin = appParams["emin"];
emax = appParams["emax"];
fovradmax = appParams["fovradmax"];
/// fovradmin = appParams["fovradmin"];
albrad = appParams["albrad"];
lonpole = appParams["lonpole"];
phasecode = appParams["phasecode"];
filtercode = appParams["filtercode"];

return true;
}

int CtsGenCommonParams::write_fits_header(fitsfile* mapFits, ProjType projection, int& status)
{
fits_update_key(mapFits, TDOUBLE, "CRVAL1", &la, NULL, &status);
fits_update_key(mapFits, TDOUBLE, "CRVAL2", &ba, NULL, &status);
char projstr1[FLEN_FILENAME];
char projstr2[FLEN_FILENAME];
if (projection == AIT) {
	strcpy(projstr1,"GLON-AIT");
	strcpy(projstr2,"GLAT-AIT");
	} else {
	strcpy(projstr1,"GLON-ARC");
	strcpy(projstr2,"GLAT-ARC");
	}
fits_update_key(mapFits, TSTRING, "CTYPE1", projstr1, NULL, &status);
fits_update_key(mapFits, TSTRING, "CTYPE2", projstr2, NULL, &status);

double xx = mdim/mres/2+0.5;
fits_update_key(mapFits, TDOUBLE, "CRPIX1", &xx, NULL, &status);
fits_update_key(mapFits, TDOUBLE, "CRPIX2", &xx, NULL, &status);	
xx = -mres;
fits_update_key(mapFits, TDOUBLE, "CDELT1", &xx, NULL, &status);
xx = mres;
fits_update_key(mapFits, TDOUBLE, "CDELT2", &xx, NULL, &status);

char unit[] = "deg";
fits_update_key(mapFits, TSTRING, "CUNIT1", unit, NULL, &status);
fits_update_key(mapFits, TSTRING, "CUNIT2", unit, NULL, &status);

char str3[] = "FK5";
fits_update_key(mapFits, TSTRING,  "RADESYS", str3, NULL, &status);

xx = 2000.0;
fits_update_key(mapFits, TDOUBLE,  "EQUINOX", &xx, NULL, &status);

fits_update_key(mapFits, TDOUBLE,  "LONPOLE", &lonpole, NULL, &status);	
WriteTime(mapFits, tmin, tmax);

fits_update_key(mapFits, TDOUBLE,  "MINENG", &emin, NULL, &status);
fits_update_key(mapFits, TDOUBLE,  "MAXENG", &emax, NULL, &status);
char str1[] = "AGILE";
fits_update_key(mapFits, TSTRING,  "TELESCOP", str1, NULL, &status);

char str2[] = "GRID";
fits_update_key(mapFits, TSTRING,  "INSTRUME", str2, NULL, &status);	

char str6[] = "T";
fits_update_key(mapFits, TSTRING,  "PIXCENT", str6, NULL, &status);	

char str7[FLEN_FILENAME] = "";
fits_update_key(mapFits, TSTRING,  "BUNIT", str7, NULL, &status);

fits_update_key(mapFits, TDOUBLE,  "FOV", &fovradmax, "Maximum off-axis angle (deg)", &status);

double fovMin;
if (GetValue("fovradmin", fovMin))
	fits_update_key(mapFits, TDOUBLE,  "FOVMIN", &fovMin, "Minimum off-axis angle (deg)", &status);
fits_update_key(mapFits, TDOUBLE,  "ALBEDO", &albrad, "Earth zenith angle (deg)", &status);
fits_update_key(mapFits, TINT,  "PHASECOD", &phasecode, "Orbital phase code", &status);

fits_update_key(mapFits, TINT,  "FILTERCO", &filtercode, "Event filter code", &status);
return status;
}

ThetaGenParams::ThetaGenParams(): CtsGenCommonParams(c_paramsTheta) {}

bool ThetaGenParams::Load(int argC, char* argV[])
{
if (!PilParams::Load(argC, argV))
	return false;
if (!AssignCommonParams())
	return false;
PilParams& appParams(*this);
const char* projstr = appParams["projection"];
if (!LoadProjection(projstr, projection))
	return false;
Print(c_paramsTheta);
return true;
}
