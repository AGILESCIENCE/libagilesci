


#include <cstring>
#include <iostream>

#include "wcstrig.h"
#include "sph.h"

#include "AgileMap.h"
#include "FitsUtils.h"

using namespace std;



void AgileMap::SetStartDate(const char* dateObs)
{
strncpy(m_dateObs, dateObs, 32);
m_dateObs[31] = 0;
}

void AgileMap::SetEndDate(const char* dateEnd)
{
strncpy(m_dateEnd, dateEnd, 32);
m_dateEnd[31] = 0;
}


void AgileMap::Copy(const AgileMap& another)
{
MatD::operator=(another);
m_lMat = another.m_lMat;
m_bMat = another.m_bMat;

m_la2 = another.m_la2;
m_ba2 = another.m_ba2;
m_lp = another.m_lp;
m_bp = another.m_bp;
m_gp = another.m_gp;
m_lonpole = another.m_lonpole;
m_emin = another.m_emin;
m_emax = another.m_emax;
m_mapIndex = another.m_mapIndex;
m_xbin = another.m_xbin;
m_x0 = another.m_x0;
m_ybin = another.m_ybin;
m_y0 = another.m_y0;

m_fovMin = another.m_fovMin;
m_fovMax = another.m_fovMax;
m_albedo = another.m_albedo;
m_phaseCode = another.m_phaseCode;
m_step = another.m_step;
strncpy(m_dateObs, another.m_dateObs, 32);
m_dateObs[31] = 0;
strncpy(m_dateEnd, another.m_dateEnd, 32);
m_dateEnd[31] = 0;
m_tstart = another.m_tstart;
m_tstop = another.m_tstop;
strncpy(m_fileName, another.m_fileName, 256);
m_fileName[255] = 0;
}



bool AgileMap::AlignedTo(const AgileMap &another) const
{
if (m_la2 == another.m_la2 && m_ba2 == another.m_ba2
		&& m_lonpole == another.m_lonpole
		&& m_xbin == another.m_xbin && m_ybin == another.m_ybin
		&& m_x0 == another.m_x0 && m_y0 == another.m_y0
		&& Cols() == another.Cols()
		&& Rows() == another.Rows())
	return true;
if (m_la2 != another.m_la2 || m_ba2 != another.m_ba2)
	cerr << "CRVAL mismatch" << endl;
if (m_lonpole != another.m_lonpole)
	cerr << "LONPOLE mismatch: " << m_lonpole << "!=" << another.m_lonpole << endl;
if (m_xbin != another.m_xbin || m_ybin != another.m_ybin)
	cerr << "CDELT mismatch" << endl;
if (m_x0 != another.m_x0 || m_y0 != another.m_y0)
	cerr << "CRPIX mismatch" << endl;
if (Cols() != another.Cols() || Rows() != another.Rows())
	cerr << "Size mismatch" << endl;
return false;
}

bool AgileMap::SameParamsAs(const AgileMap &another) const
{
if (AlignedTo(another))
	if (m_emin==another.m_emin && m_emax==another.m_emax)
		return true;
if (m_emin != another.m_emin || m_emax != another.m_emax)
	cerr << "Energy mismatch" << endl;
return false;
}


int AgileMap::Read(const char* fileName)
{
m_fileName[0] = 0;

FitsFile f;
if (!f.Open(fileName)) {
	cerr << "ERROR " << f.Status() << " opening " << fileName << endl;
	return f.Status();
	}

strncpy(m_fileName, fileName, 1024);
m_fileName[1024] = 0;
int bitpix, naxis;
long naxes[2] = { 1, 1 };
if (fits_get_img_param(f, 2, &bitpix, &naxis, naxes, f)) {
	cerr << "Wrong image parameters" << endl;
	return 1;
	}
if (naxis != 2) {
	cerr << "Error: only 2D images are supported" << endl;
	return 1;
	}

/// zzz MatD mat(2, naxes[1], naxes[0]);
MatD mat(naxes[1], naxes[0]);
long fpixel[2] = { 1, 1 };
fits_read_pix(f, TDOUBLE, fpixel, mat.Size(), NULL, mat.Buffer(), NULL, f);
mat.TransposeTo(*this);

/// Mandatory keywords
f.ReadKey("CDELT1", &m_xbin);
f.ReadKey("CRPIX1", &m_x0);
f.ReadKey("CDELT2", &m_ybin);
f.ReadKey("CRPIX2", &m_y0);
f.ReadKey("CRVAL1", &m_la2);
f.ReadKey("CRVAL2", &m_ba2);
f.ReadKey("LONPOLE", &m_lonpole);
f.ReadKey("MINENG", &m_emin);
/// Optional keywords
f.ReadKey("MAXENG", &m_emax, 99999);
f.ReadKey("INDEX", &m_mapIndex, -2.1);
f.ReadKey("SC-Z-LII", &m_lp, double(-999));
f.ReadKey("SC-Z-BII", &m_bp, double(-999));
f.ReadKey("SC-LONPL", &m_gp, 0.0);
f.ReadKey("DATE-OBS", m_dateObs, "");
f.ReadKey("DATE-END", m_dateEnd, "");
f.ReadKey("TSTART", &m_tstart, 0.0);
f.ReadKey("TSTOP", &m_tstop, 0.0);
f.ReadKey("FOVMIN", &m_fovMin, 0.0);
f.ReadKey("FOV", &m_fovMax, 0.0);
f.ReadKey("ALBEDO", &m_albedo, 0.0);
f.ReadKey("PHASECOD", &m_phaseCode, 0l);
f.ReadKey("STEP", &m_step, 0.0);

if (f.Status())
	cerr << "ERROR " << f.Status() << " reading " << fileName << endl;
else
	Eval_lb();
return f.Status();
}



int AgileMap::Write(const char* fileName) const
{
FitsFile f;
if (!f.Create(fileName)) {
	cerr << "ERROR " << f.Status() << " creating " << fileName << endl;
	return f.Status();
	}

MatD mat;
TransposeTo(mat);

int bitpix = DOUBLE_IMG;
int naxis = 2;
long naxes[2] = { mat.Cols(), mat.Rows() };
long fpixel[2] = { 1, 1 };
fits_create_img(f, bitpix, naxis, naxes, f);
fits_write_pix(f, TDOUBLE, fpixel, mat.Size(), const_cast<double*>(mat.Buffer()), f);

f.WriteKey("CTYPE1", "GLON-ARC");
f.WriteKey("CTYPE2", "GLAT-ARC");
f.WriteKey("CRPIX1", m_x0);
f.WriteKey("CRVAL1", m_la2);
f.WriteKey("CDELT1", m_xbin);
f.WriteKey("CRPIX2", m_y0);
f.WriteKey("CRVAL2", m_ba2);
f.WriteKey("CDELT2", m_ybin);
f.WriteKey("LONPOLE", m_lonpole);
f.WriteKey("MINENG", m_emin);
f.WriteKey("MAXENG", m_emax);
f.WriteKey("INDEX", m_mapIndex);
f.WriteKey("SC-Z-LII", m_lp);
f.WriteKey("SC-Z-BII", m_bp);
f.WriteKey("SC-LONPL", m_gp);

if (m_dateObs[0])
	f.WriteKey("DATE-OBS", m_dateObs);
if (m_dateEnd[0])
	f.WriteKey("DATE-END", m_dateEnd);

f.WriteKey("TSTART", m_tstart);
f.WriteKey("TSTOP", m_tstop);
f.WriteKey("FOVMIN", m_fovMin);
f.WriteKey("FOV", m_fovMax);
f.WriteKey("ALBEDO", m_albedo);
f.WriteKey("PHASECOD", m_phaseCode);

if (m_step)
	f.WriteKey("STEP", m_step);

if (f.Status())
	cerr << "ERROR " << f.Status() << " writing to " << fileName << endl;
return f.Status();
}

int AgileMap::WriteWithAllMetadata(const char* fileName) const
{
FitsFile f;
if (!f.Create(fileName)) {
	cerr << "ERROR " << f.Status() << " creating " << fileName << endl;
	return f.Status();
	}

MatD mat;
TransposeTo(mat);

FitsFile from;
if (!from.Open(m_fileName)) {
    cerr << "ERROR " << f.Status() << " reopening " << m_fileName << endl;
    return f.Status();
}

int status = 0;
fits_copy_header(from, f, &status);
if (status) {
    cerr << "ERROR " << status << " copying hdu from " << m_fileName << " to " << fileName << endl;
    return status;
}
from.Close();

f.UpdateKey("NAXIS", 2);
f.UpdateKey("NAXIS1", mat.Cols());
f.UpdateKey("NAXIS2", mat.Rows());

long fpixel[2] = { 1, 1 };
fits_write_pix(f, TDOUBLE, fpixel, mat.Size(), const_cast<double*>(mat.Buffer()), f);

f.UpdateKey("CTYPE1", "GLON-ARC");
f.UpdateKey("CTYPE2", "GLAT-ARC");
f.UpdateKey("CRPIX1", m_x0);
f.UpdateKey("CRVAL1", m_la2);
f.UpdateKey("CDELT1", m_xbin);
f.UpdateKey("CRPIX2", m_y0);
f.UpdateKey("CRVAL2", m_ba2);
f.UpdateKey("CDELT2", m_ybin);
f.UpdateKey("LONPOLE", m_lonpole);
f.UpdateKey("MINENG", m_emin);
f.UpdateKey("MAXENG", m_emax);
f.UpdateKey("INDEX", m_mapIndex);
f.UpdateKey("SC-Z-LII", m_lp);
f.UpdateKey("SC-Z-BII", m_bp);
f.UpdateKey("SC-LONPL", m_gp);

if (m_dateObs[0])
	f.UpdateKey("DATE-OBS", m_dateObs);
if (m_dateEnd[0])
	f.UpdateKey("DATE-END", m_dateEnd);

f.UpdateKey("TSTART", m_tstart);
f.UpdateKey("TSTOP", m_tstop);
f.UpdateKey("FOVMIN", m_fovMin);
f.UpdateKey("FOV", m_fovMax);
f.UpdateKey("ALBEDO", m_albedo);
f.UpdateKey("PHASECOD", m_phaseCode);

if (m_step)
	f.UpdateKey("STEP", m_step);

if (f.Status())
	cerr << "ERROR " << f.Status() << " writing to " << fileName << endl;
return f.Status();
}



double AgileMap::Area(int i, int j) const
{
double pixel = DEG2RAD * DEG2RAD * fabs(m_xbin * m_ybin);
return pixel * Sinaa(DEG2RAD*theta(i,j));
}

double AgileMap::phi(int i, int j) const
{
// return atan2d(x(i,j),-y(i,j));
return atan2d(x(i),-y(j));
}


bool AgileMap::GetRowCol(double lDeg, double bDeg, int* row, int* col) const
{
double la = GetMapCenterL() * DEG2RAD;
double ba = GetMapCenterB() * DEG2RAD;
double cosba = cos(ba);
double sinba = sin(ba);

double mres = fabs(GetXbin());
/// double half_mdim = (GetX0() - 0.5)*mres;
double half_mapX = Cols()/2;
double half_mapY = Rows()/2;

/// Minimize round-off errors
while (lDeg<0)
	lDeg += 360;
while (lDeg>=360)
	lDeg -= 360;

double l = lDeg*DEG2RAD;
double b = bDeg*DEG2RAD;
double theta = sin(b)*sinba+cos(b)*cosba*cos(l-la);

double x, y;
if(theta == 1.0) {
    x = 0.0;
    y = 0.0;
}
else {
    if (theta < -1.0)
        theta = M_PI;
    else if (theta > 1.0)
        theta = 0.0;
    else
        theta = acos(theta);
    x = RAD2DEG/Sinaa(theta) * cos(b)*sin(l-la);
    y = RAD2DEG/Sinaa(theta) * (sin(b)*cosba - cos(b)*sinba*cos(l-la));
}

*row = int(floor(half_mapX-x/mres));
*col = int(floor(half_mapY+y/mres));
return *row>=0 && *row<Rows() && *col>=0 && *col<Cols();
}


void AgileMap::GetCoords(int row, int col, double* lDeg, double* bDeg) const
{
double eul[5];
eul[0] = m_la2;
eul[1] = 90.0 - m_ba2;
eul[2] = m_lonpole;
eul[3] = cosd(eul[1]);
eul[4] = sind(eul[1]);
double x = (row+1-m_x0)*m_xbin;
double y = (col+1-m_y0)*m_ybin;
double theta = 90.0 - sqrt(x*x+y*y);
double phi = atan2d(x,-y);
sphx2s(eul, 1, 1, 0, 0, &phi, &theta, lDeg, bDeg);
}


MatD AgileMap::SrcDist(double lng, double lat) const
{
MatD xxx(Rows(), Cols());
for (int i=0; i<Rows(); ++i)
	for (int j=0; j< Cols(); ++j)
		xxx(i,j) = SrcDist(i, j, lng, lat);
return xxx;
}



void AgileMap::Eval_lb()
{
double eul[5];
eul[0] = m_la2;
eul[1] = 90.0 - m_ba2;
eul[2] = m_lonpole;
eul[3] = cosd(eul[1]);
eul[4] = sind(eul[1]);
int rows = Rows();
int cols = Cols();

/** zzz was ReshapeTo */
m_lMat.ResizeTo(*this);
m_bMat.ResizeTo(*this);

double l, b;
for (int i=0; i<rows; ++i)
	for (int j=0; j<cols; ++j) {
		double theta2 = 90.0 - theta(i,j);
		double phi2 = phi(i,j);
		sphx2s(eul, 1, 1, 0, 0, &phi2, &theta2, &l, &b);
		m_lMat(i,j) = l;
		m_bMat(i,j) = b;
		}
}

