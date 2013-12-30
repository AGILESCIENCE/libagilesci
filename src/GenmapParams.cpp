

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
char * secstr = "s";

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




bool Maplist::Store(const char* maplist)
{
	bool stillgoing = true;
	bool succeeded = false;
	MaplistEntry mapspec;
	fstream mapstream(maplist);
	while (stillgoing) {
		mapstream >> mapspec.fovradmin >> mapspec.fovradmax >> mapspec.emin >> mapspec.emax >> mapspec.index;
		if (mapstream.good()) {
			push_back(mapspec);
			succeeded = true;
		} else
			stillgoing = false;
	}
	return succeeded;
}

void Maplist::Print() // const
{
cout << "Energy ranges (FOVmin [deg], FOVmax [deg], Emin [MeV], Emax [MeV], spectral index) :" << endl;
/// for (const_iterator i = cbegin(); i != cend(); i++)
for (iterator i = begin(); i != end(); i++)
	cout << setprecision(2) << "     "
		<< i->fovradmin << "   "
		<< i->fovradmax << "     "
		<< i->emin << "   "
		<< i->emax << "   "
		<< i->index << endl;
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



const PilDescription c_paramsExp[] = {
	{ PilString, "outfile", "Output file name" },
	{ PilString, "logfile", "Grid log index file name" },
	{ PilString, "sarFileName", "Effective area file name" },
	{ PilString, "edpFileName", "Energy dispersion file name" },
	{ PilString, "maplist", "Maplist file name" },
	{ PilString, "timelist", "Time intervals file name" },
	{ PilReal, "mdim", "Size of Map (degrees)" },
	{ PilReal, "mres", "Bin size (degrees)" },
	{ PilReal, "la", "Longitude of map center (galactic)" },
	{ PilReal, "ba", "Latitude of map center (galactic)" },
	{ PilReal, "lonpole", "Rotation of map (degrees)" },
	{ PilReal, "albrad", "Radius of earth albedo (degrees)" },
	{ PilReal, "y_tol", "Boresight movement tolerance (degrees)" },
	{ PilReal, "roll_tol", "Roll tolerance (degrees)" },
	{ PilReal, "earth_tol", "Earth tolerance (degrees)" },
	{ PilInt, "phasecode", "Orbital phase code" },
	{ PilString, "projection", "Projection (ARC or AIT)" },
	{ PilReal, "binstep", "Bin step size" },
	{ PilReal, "timestep", "LOG file step size" },
	{ PilReal, "index", "Spectral index" },
	{ PilReal, "tmin", "Initial time (sec)" },
	{ PilReal, "tmax", "Final time (sec)" },
	{ PilReal, "emin", "Minimum energy" },
	{ PilReal, "emax", "Maximum energy" },
	{ PilReal, "fovradmin", "Min radius of field of view (degrees)" },
	{ PilReal, "fovradmax", "Max radius of field of view (degrees)" },
	{ PilNone, "", "" }
	};


/// Same as above except:
/// Removed: edpFileName maplist mres binstep
/// Added: timeslot
const PilDescription c_paramsExpT[] = {
	{ PilString, "outfile", "Output file name" },
	{ PilString, "logfile", "Grid log index file name" },
	{ PilString, "sarFileName", "Effective area file name" },
	{ PilString, "timelist", "Time intervals file name" },
	{ PilReal, "mdim", "Size of Map (degrees)" },
	{ PilReal, "la", "Longitude of map center (galactic)" },
	{ PilReal, "ba", "Latitude of map center (galactic)" },
	{ PilReal, "lonpole", "Rotation of map(degrees)" },
	{ PilReal, "albrad", "Radius of earth albedo (degrees)" },
	{ PilReal, "y_tol", "Boresight movement tolerance (degrees)" },
	{ PilReal, "roll_tol", "Roll tolerance (degrees)" },
	{ PilReal, "earth_tol", "Earth tolerance (degrees)" },
	{ PilInt, "phasecode", "Orbital phase code" },
	{ PilString, "projection", "Projection (ARC or AIT)" },
	{ PilReal, "timestep", "LOG file step size" },
	{ PilReal, "timeslot", "Time slot" },
	{ PilReal, "index", "Spectral index" },
	{ PilReal, "tmin", "Initial time (sec)" },
	{ PilReal, "tmax", "Final time (sec)" },
	{ PilReal, "emin", "Minimum energy" },
	{ PilReal, "emax", "Maximum energy" },
	{ PilReal, "fovradmin", "Min off-axis angle (degrees)" },
	{ PilReal, "fovradmax", "Max off-axis angle (degrees)" },
	{ PilNone, "", "" }
	};



bool ExpGenCommonParams::AssignCommonParams()
{
PilParams& appParams(*this);
outfile = appParams["outfile"];
logfile = appParams["logfile"];
sarFileName = appParams["sarFileName"];
timelist = appParams["timelist"];
mdim = appParams["mdim"];
la = appParams["la"];
ba = appParams["ba"];
lonpole = appParams["lonpole"];
albrad = appParams["albrad"];
y_tol = appParams["y_tol"];
roll_tol = appParams["roll_tol"];
earth_tol = appParams["earth_tol"];
phasecode = appParams["phasecode"];

const char* projstr = appParams["projection"];
if (!LoadProjection(projstr, projection))
	return false;

timestep = appParams["timestep"];
index = appParams["index"];
tmin = appParams["tmin"];
tmax = appParams["tmax"];
emin = appParams["emin"];
emax = appParams["emax"];
fovradmin = appParams["fovradmin"];
fovradmax = appParams["fovradmax"];

if (!LoadTimeList(timelist, intervals, tmin, tmax))
	return false;
return true;
}


ExpGenParams::ExpGenParams(): ExpGenCommonParams(c_paramsExp) {}

bool ExpGenParams::Load(int argC, char* argV[])
{
if (!PilParams::Load(argC, argV))
	return false;

if (!AssignCommonParams())
	return false;
PilParams& appParams(*this);
edpFileName = appParams["edpFileName"];
maplist = appParams["maplist"];
mres = appParams["mres"];
binstep = appParams["binstep"];
mxdim = long(mdim/mres+0.1);

if (strcmp(maplist, "None")) {
	maps.Store(maplist);
	maps.Print();
	}
else {
	MaplistEntry mapspec;
	mapspec.fovradmin = fovradmin;
	mapspec.fovradmax = fovradmax;
	mapspec.emin = emin;
	mapspec.emax = emax;
	mapspec.index= index;
	maps.push_back(mapspec);
	}

PilParams::Print(c_paramsExp);
return true;
}


ExpGenParamsT::ExpGenParamsT(): ExpGenCommonParams(c_paramsExpT) {}

bool ExpGenParamsT::Load(int argC, char* argV[])
{
if (!PilParams::Load(argC, argV))
	return false;

if (!AssignCommonParams())
	return false;
PilParams& appParams(*this);
timeslot = appParams["timeslot"];

PilParams::Print(c_paramsExpT);
return true;
}







static const char* AfterSlashPtr(const char* str)
{
const char* slash = strrchr(str, '/');
return slash ? slash+1 : str;
}

int ExpGenParams::write_fits_header(long i, fitsfile *mapFits, int & status)
{
	
//	std::cout<< "writing header........................................" << std::endl<< std::endl;	

	double fovradmin = maps[i].fovradmin;
	double fovradmax = maps[i].fovradmax;
	double emin = maps[i].emin;
	double emax = maps[i].emax;
	double spectral_index = maps[i].index;
	
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

	char strValue[72];
	strValue[68] = 0;
	strncpy(strValue, AfterSlashPtr(sarFileName), 68);
	fits_update_key(mapFits, TSTRING,  "SARFILE", strValue, NULL, &status);
	strncpy(strValue, AfterSlashPtr(edpFileName), 68);
	fits_update_key(mapFits, TSTRING,  "EDPFILE", strValue, NULL, &status);

	char str3[] =  "FK5";
	fits_update_key(mapFits, TSTRING,  "RADESYS", str3, NULL, &status);
	
	xx = 2000.;
	fits_update_key(mapFits, TDOUBLE,  "EQUINOX", &xx, NULL, &status);

	fits_update_key(mapFits, TDOUBLE,  "LONPOLE", &lonpole, NULL, &status);
	WriteTime(mapFits, tmin, tmax);
	
	fits_update_key(mapFits, TDOUBLE,  "MINENG", &emin, NULL, &status);
	fits_update_key(mapFits, TDOUBLE,  "MAXENG", &emax, NULL, &status);
	char * str1 =  "AGILE";		
	fits_update_key(mapFits, TSTRING,  "TELESCOP", str1, NULL, &status);

	char * str2 =  "GRID";
	fits_update_key(mapFits, TSTRING,  "INSTRUME", str2, NULL, &status);	
	
	char * str6 =  "T";
	
	fits_update_key(mapFits, TSTRING,  "PIXCENT", str6, NULL, &status);	

	char str7[FLEN_FILENAME];
		strcpy(str7, "cm**2 s sr");
		fits_update_key(mapFits, TDOUBLE,  "YTOL", &y_tol, "Pointing direction tolerance (deg)", &status);
		fits_update_key(mapFits, TDOUBLE,  "ROLTOL", &roll_tol, "Roll angle tolerance (deg)", &status);
		fits_update_key(mapFits, TDOUBLE,  "EARTOL", &earth_tol, "Pointing direction tolerance (deg)", &status);
		fits_update_key(mapFits, TINT,  "STEP", &binstep, "Map interpolation step size", &status);
		fits_update_key(mapFits, TINT,  "TIMESTEP", &timestep, "Log file step size", &status);

	fits_update_key(mapFits, TSTRING,  "BUNIT", str7, NULL, &status);
	
	fits_update_key(mapFits, TDOUBLE,  "FOV", &fovradmax, "Radius of field of view (deg)", &status);
	fits_update_key(mapFits, TDOUBLE,  "FOVMIN", &fovradmin, "Minimum off-axis angle (deg)", &status);
	fits_update_key(mapFits, TDOUBLE,  "ALBEDO", &albrad, "Earth zenith angle (deg)", &status);
	fits_update_key(mapFits, TINT,  "PHASECOD", &phasecode, "Orbital phase code", &status);
	
		fits_update_key(mapFits, TDOUBLE, "INDEX", &spectral_index, NULL, &status);

	return status;
}





const PilDescription c_paramsCts[] = {
	{ PilString, "outfile", "Output file name" },
	{ PilString, "evtfile", "Event file index file name" },
	{ PilString, "timelist", "Time intervals list" },
	{ PilReal, "mdim", "Size of Map (degrees)" },
	{ PilReal, "mres", "Bin size (degrees)" },
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
	{ PilReal, "fovradmin", "Min off-axis angle (degrees)" },
	{ PilReal, "fovradmax", "Max off-axis angle (degrees)" },
	{ PilNone, "", "" }
	};


const PilDescription c_paramsCtsT[] = {
	{ PilString, "outfile", "Output file name" },
	{ PilString, "evtfile", "Event file index file name" },
	{ PilString, "timelist", "Time intervals list" },
	{ PilReal, "mdim", "Selection diameter (degrees)" },
	{ PilReal, "mres", "Bin size (degrees)" },
	{ PilReal, "la", "Longitude of selection center (galactic)" },
	{ PilReal, "ba", "Latitude of selection center (galactic)" },
	{ PilReal, "lonpole", "Rotation of map(degrees)" },
	{ PilReal, "albrad", "Radius of earth albedo (degrees)" },
	{ PilInt, "phasecode", "Orbital phase code" },
	{ PilInt, "filtercode", "Event filter code" },
	{ PilReal, "tmin", "Initial time(sec)" },
	{ PilReal, "tmax", "Final time(sec)" },
	{ PilReal, "emin", "Min energy" },
	{ PilReal, "emax", "Max energy" },
	{ PilReal, "fovradmin", "Min off-axis angle (degrees)" },
	{ PilReal, "fovradmax", "Max off-axis angle (degrees)" },
	{ PilReal, "timeslot", "Time slot" },
	{ PilNone, "", "" }
	};

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





CtsGenParams::CtsGenParams(): CtsGenCommonParams(c_paramsCts) {}


bool CtsGenParams::Load(int argC, char* argV[])
{
if (!PilParams::Load(argC, argV))
	return false;
if (!AssignCommonParams())
	return false;
PilParams& appParams(*this);
timelist = appParams["timelist"];
fovradmin = appParams["fovradmin"];
if (!LoadTimeList(timelist, intervals, tmin, tmax))
	return false;
const char* projstr = appParams["projection"];
if (!LoadProjection(projstr, projection))
	return false;
Print(c_paramsCts);
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

char * unit =  "deg";
fits_update_key(mapFits, TSTRING, "CUNIT1", unit, NULL, &status);
fits_update_key(mapFits, TSTRING, "CUNIT2", unit, NULL, &status);

char * str3 =  "FK5";
fits_update_key(mapFits, TSTRING,  "RADESYS", str3, NULL, &status);

xx = 2000.0;
fits_update_key(mapFits, TDOUBLE,  "EQUINOX", &xx, NULL, &status);

fits_update_key(mapFits, TDOUBLE,  "LONPOLE", &lonpole, NULL, &status);	
WriteTime(mapFits, tmin, tmax);

fits_update_key(mapFits, TDOUBLE,  "MINENG", &emin, NULL, &status);
fits_update_key(mapFits, TDOUBLE,  "MAXENG", &emax, NULL, &status);
char * str1 =  "AGILE";
fits_update_key(mapFits, TSTRING,  "TELESCOP", str1, NULL, &status);

char * str2 =  "GRID";
fits_update_key(mapFits, TSTRING,  "INSTRUME", str2, NULL, &status);	

char * str6 =  "T";
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




CtsGenParamsT::CtsGenParamsT(): CtsGenCommonParams(c_paramsCtsT) {}


bool CtsGenParamsT::Load(int argC, char* argV[])
{
if (!PilParams::Load(argC, argV))
	return false;
if (!AssignCommonParams())
	return false;
PilParams& appParams(*this);
timelist = appParams["timelist"];
fovradmin = appParams["fovradmin"];
timeslot = appParams["timeslot"];
Print(c_paramsCtsT);
return true;
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













#ifdef _DEAD_CODE_
string ExpGenParams::evtexpr(long i)
{
	char expr[1024];
	sprintf(expr, "TIME >= %f && TIME <= %f && ENERGY >= %f && ENERGY <= %f && PH_EARTH > %f && THETA < %f && THETA >= %f", 
		  tmin, tmax, maps[i].emin, maps[i].emax, albrad, maps[i].fovradmax, maps[i].fovradmin);
	if ((phasecode & 1) == 1) strcat(expr, " && PHASE .NE. 0");
	if ((phasecode & 2) == 2) strcat(expr, " && PHASE .NE. 1");
	if ((phasecode & 4) == 4) strcat(expr, " && PHASE .NE. 2");
	if ((phasecode & 8) == 8) strcat(expr, " && PHASE .NE. 3");
	if ((phasecode & 16) == 16) strcat(expr, " && PHASE .NE. 4");
	if ((filtercode & 1) == 1) strcat(expr, " && EVSTATUS .NE. 'L'");
	if ((filtercode & 2) == 2) strcat(expr, " && EVSTATUS .NE. 'G'");
	if ((filtercode & 4) == 4) strcat(expr, " && EVSTATUS .NE. 'S'");
/*	if (filtercode == 1)
		strcat(expr, " && EVSTATUS == 'G'");
	else if (filtercode == 2)
		strcat(expr, " && EVSTATUS == 'L'");*/
	return string(expr);
}

string ExpGenParams::logexpr()
{
	char expr[1024];
	sprintf(expr, "TIME >= %f && TIME <= %f && LIVETIME > 0 && LOG_STATUS == 0 && MODE == 2", tmin, tmax);
	if ((phasecode & 1) == 1) strcat(expr, " && PHASE .NE. 0");
	if ((phasecode & 2) == 2) strcat(expr, " && PHASE .NE. 1");
	if ((phasecode & 4) == 4) strcat(expr, " && PHASE .NE. 2");
	if ((phasecode & 8) == 8) strcat(expr, " && PHASE .NE. 3");
	if ((phasecode & 16) == 16) strcat(expr, " && PHASE .NE. 4");
	return string(expr);
}

#endif


