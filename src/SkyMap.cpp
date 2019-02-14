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

#include <cstring>
#include <ctime>
#include <iostream>
#include <fitsio.h>


#include "SkyMap.h"

using std::cerr;
using std::cout;
using std::endl;


static bool ReadFitsKey(fitsfile* fptr, const char* name, float* value, float defaultVal, int* status)
{
if (!*status) {
	int status2 = 0;
	fits_read_key(fptr, TFLOAT, const_cast<char*>(name), value, NULL, &status2);
	if (status2==0)
		return true;
	else if (status2=KEY_NO_EXIST)
		*value = defaultVal;
	else
		*status = status2;
	}
return false;
}

static void ReadFitsKey(fitsfile* fptr, const char* name, float* value, int* status)
{
fits_read_key(fptr, TFLOAT, const_cast<char*>(name), value, NULL, status);
}


static void ReadFitsKey(fitsfile* fptr, const char* name, char* value, int* status)
{
fits_read_key(fptr, TSTRING, const_cast<char*>(name), value, NULL, status);
}


static void CommentKeyword(const char* kwName, int status)
{
if (status)
	cerr << "Problem writing " << kwName << std::endl;
}

static void WriteKeyword(fitsfile* fptr, int* status, const char* kwName, bool value, const char* comment = 0)
{
fits_write_key(fptr, TLOGICAL, const_cast<char*>(kwName), &value, const_cast<char*>(comment), status);
CommentKeyword(kwName, *status);
}

static void WriteKeyword(fitsfile* fptr, int* status, const char* kwName, int value, const char* comment = 0)
{
fits_write_key(fptr, TINT, const_cast<char*>(kwName), &value, const_cast<char*>(comment), status);
CommentKeyword(kwName, *status);
}

static void WriteKeyword(fitsfile* fptr, int* status, const char* kwName, float value, const char* comment = 0)
{
fits_write_key(fptr, TFLOAT, const_cast<char*>(kwName), &value, const_cast<char*>(comment), status);
CommentKeyword(kwName, *status);
}

static void WriteKeyword(fitsfile* fptr, int* status, const char* kwName, const char* value, const char* comment = 0)
{
fits_write_key(fptr, TSTRING, const_cast<char*>(kwName), const_cast<char*>(value), const_cast<char*>(comment), status);
CommentKeyword(kwName, *status);
}





void SkyMap::Print() const
{
cout << "m_rows = " << m_rows << endl;
cout << "m_cols = " << m_cols << endl;
cout << "m_xRes = " << m_xRes << endl;
cout << "m_yRes = " << m_yRes << endl;
cout << "m_xRef = " << m_xRef << endl;
cout << "m_yRef = " << m_yRef << endl;
cout << "m_xPoint = " << m_xPoint << endl;
cout << "m_yPoint = " << m_yPoint << endl;
cout << "m_minX = " << m_minX << endl;
cout << "m_maxX = " << m_maxX << endl;
cout << "m_minY = " << m_minY << endl;
cout << "m_maxY = " << m_maxY << endl;
for (int i=0; i<m_kwrdCount; ++i)
	cout << m_kwrdText+m_kwrdIndices[i] << endl;
}



void SkyMap::SetKeywords(const char* const keywords[])
{
delete[] m_kwrdIndices;
delete[] m_kwrdText;
m_kwrdIndices = 0;
m_kwrdText = 0;
m_kwrdCount = 0;
int size = 0;
while (keywords && keywords[m_kwrdCount])
	size += strlen(keywords[m_kwrdCount++]);
if (!m_kwrdCount)
	return;
m_kwrdIndices = new int[m_kwrdCount];
m_kwrdText = new char[size+m_kwrdCount];
int offset = 0;
for (int i=0; i<m_kwrdCount; ++i) {
	m_kwrdIndices[i] = offset;
	strcpy(m_kwrdText+offset, keywords[i]);
	offset += 1+strlen(keywords[i]);
	}
}



void SkyMap::CopyKeywords(const SkyMap& other)
{
delete[] m_kwrdIndices;
delete[] m_kwrdText;
m_kwrdIndices = 0;
m_kwrdText = 0;
m_kwrdCount = other.m_kwrdCount;
if (m_kwrdCount) {
	int lastIndex = other.m_kwrdIndices[m_kwrdCount-1];
	int otherSize = lastIndex+strlen(other.m_kwrdText+lastIndex)+1;
	m_kwrdIndices = new int[m_kwrdCount];
	m_kwrdText = new char[otherSize];
	memcpy(m_kwrdIndices, other.m_kwrdIndices, sizeof(int)*m_kwrdCount);
	memcpy(m_kwrdText, other.m_kwrdText, otherSize);
	}
}


void SkyMap::Copy(const SkyMap& other)
{
delete[] m_values;
m_rows = other.m_rows;
m_cols = other.m_cols;
int count = m_rows*m_cols;
if (count) {
	m_values = new float[count];
	memcpy(m_values, other.m_values, sizeof(float)*count);
	}
else
	m_values = 0;
m_xRes = other.m_xRes;
m_yRes = other.m_yRes;
m_xRef = other.m_xRef;
m_yRef = other.m_yRef;
m_xPoint = other.m_xPoint;
m_yPoint = other.m_yPoint;

m_minX = other.m_minX;
m_maxX = other.m_maxX;
m_minY = other.m_minY;
m_maxY = other.m_maxY;

CopyKeywords(other);
}



void SkyMap::Resize(int rows, int cols)
{
m_rows = rows;
m_cols = cols;
delete[] m_values;

delete[] m_kwrdIndices;
delete[] m_kwrdText;
m_kwrdCount = 0;
m_kwrdIndices = 0;
m_kwrdText = 0;

m_values = new float[m_rows*m_cols];
Clear();
}


SkyMap::SkyMap(int rows, int cols, float xRes, float yRes, float xPoint, float yPoint, float xRef, float yRef)
{
Reset();
Resize(rows, cols);
m_xRes = xRes;
m_yRes = yRes;
//m_xRef = xRef>=0 ? xRef : float(rows+1)/2;
//m_yRef = yRef>=0 ? yRef : float(cols+1)/2;
m_xRef = xRef;
m_yRef = yRef;
m_xPoint = xPoint;
m_yPoint = yPoint;
SetLimits();
}




void SkyMap::SetAllValues(float value)
{
int count = m_rows*m_cols;
for (int i=0; i<count; ++i)
	m_values[i] = value;
}


void SkyMap::GetRefPoint(float* xRef, float* yRef, float* xPoint, float* yPoint) const
{
if (xRef)
	*xRef = m_xRef;
if (yRef)
	*yRef = m_yRef;
if (xPoint)
	*xPoint = m_xPoint;
if (yPoint)
	*yPoint = m_yPoint;
}


void SkyMap::SetRefPoint(float xRef, float yRef, float xPoint, float yPoint)
{
m_xRef = xRef;
m_yRef = yRef;
m_xPoint = xPoint;
m_yPoint = yPoint;
SetLimits();
}



namespace {

/// refBin may indicate either the centre or the side of the bin
#define REFPOINT_CENTER

	float PointOf(float refPoint, float refBin, float res, int bin)
	{
#ifdef REFPOINT_CENTER
	return refPoint+(0.5-refBin+bin)*res;
#else
	return refPoint+(1-refBin+bin)*res;
#endif
	}
	
	int BinOf(float refPoint, float refBin, float res, float point)
	{
#ifdef REFPOINT_CENTER
	return int(refBin-1+(point-refPoint)/res);
#else
	return int(refBin-0.5+(point-refPoint)/res);
#endif
	}

	void FloatSwap(float& a, float& b)
	{
	float t = a;
	a = b;
	b = t;
	}

	void EvalRange(float point, float ref, float res, int length, float& min, float& max)
	{
	/// float imin = point+(1-ref)*res;
	/// float imax = point+(1-ref+length)*res;
	min = PointOf(point, ref, res, 0);
	max = PointOf(point, ref, res, length);
	if (res<0)
		FloatSwap(min, max);
	/// cerr << "point=" << point << " ref=" << ref << " res=" << res << " length=" << length << " min=" << min << " max=" << max << endl;
	}

}




void SkyMap::SetLimits()
{
EvalRange(m_xPoint, m_xRef, m_xRes, m_cols, m_minX, m_maxX);
EvalRange(m_yPoint, m_yRef, m_yRes, m_rows, m_minY, m_maxY);
}



void SkyMap::Range(float* minX, float* maxX, float* minY, float* maxY) const
{
if (minX)
	*minX = m_minX;
if (maxX)
	*maxX = m_maxX;
if (minY)
	*minY = m_minY;
if (maxY)
	*maxY = m_maxY;
}


static bool InsideAngle(float& angle, float min, float max)
{
while (angle<min)
	angle += 360;
while (angle>max)
	angle -= 360;
return angle>=min && angle<=max;
}



bool SkyMap::Point2Bin(float x, float y, int* row, int* col) const
{
bool inside = true;
if (row) {
	float minY, maxY;
	Range(0, 0, &minY, &maxY);
	inside = InsideAngle(y, minY, maxY);
	*row = BinOf(m_yPoint, m_yRef, m_yRes, y);
	}
if (col) {
	float minX, maxX;
	Range(&minX, &maxX, 0, 0);
	inside = inside && InsideAngle(x, minX, maxX);
	*col = BinOf(m_xPoint, m_xRef, m_xRes, x);
	}
return inside;
}

bool SkyMap::Bin2Point(int row, int col, float* x, float* y) const
{
if (x)
	*x = m_xPoint+(-m_xRef+col+1)*m_xRes;
if (y)
	*y = m_yPoint+(-m_yRef+row+1)*m_yRes;
return row>=0 && row<m_rows && col>=0 && col<m_cols;
}



/**
static void PrintKeywords(fitsfile *fptr, int* status)
{
int keyCount;
fits_get_hdrspace(fptr, &keyCount, 0, status);
cout << keyCount << " keywords" << endl;
char keyname[128];
char value[128];
char comment[128];
for (int i=0; i<keyCount; ++i) {
	fits_read_keyn(fptr, i+1, keyname, value, comment, status);
	cout << i+1 << " " << keyname << " " << value << " " << comment << endl;
	}
}
*/



namespace
{
/// #define MORE_KEYWORDS

	const char* const knownKeywords[] = {
		"XTENSION",
		"BITPIX",
		"NAXIS",
		"NAXIS1",
		"NAXIS2",
		"PCOUNT",
		"GCOUNT",
		"BSCALE",
		"BZERO",
		"CDELT1",
		"CRPIX1",
		"CRVAL1",
		"CDELT2",
		"CRPIX2",
		"CRVAL2",
		"EXTNAME",
#ifdef MORE_KEYWORDS
		"CTYPE1",
		"CUNIT1",
		"CTYPE2",
		"CUNIT2",
		"PIXCENT",
		"BUNIT",
		"PRIMTYPE",
		"INSTRUME",
#endif
		0};


	bool IsUnknownKeyword(const char* keyword)
	{
	for (int i=0; knownKeywords[i]; ++i)
		if (!strcmp(keyword, knownKeywords[i]))
			return false;
	return true;
	}
	
	void ReadUnknKeywords(fitsfile *fptr, int& count, int*& offsets, char*& text, int* status)
	{
	count = 0;
	offsets = 0;
	text = 0;
	if (*status)
		return;
	
	int l_count;
	fits_get_hdrspace(fptr, &l_count, 0, status);
	if (l_count==0 || *status)
		return;
	char* l_text = new char[l_count*84];
	int* l_offsets = new int[l_count];
	int index = 0;
	int offset = 0;
	for (int i=0; i<l_count && !*status; ++i) {
		char keyname[128];
		char value[128];
		fits_read_keyn(fptr, i+1, keyname, value, 0, status);
		if (IsUnknownKeyword(keyname)) {
			char record[128];
			fits_read_record(fptr, i+1, record, status);
			int length = strlen(record);
			strcpy(l_text+offset, record);
			l_offsets[index] = offset;
			offset += length+1;
			++index;
			}
		}
	
	if (index && !*status) {
		count = index;
		offsets = new int[count];
		memcpy(offsets, l_offsets, sizeof(int)*count);
		text = new char[offset];
		memcpy(text, l_text, offset);
		}
	delete[] l_offsets;
	delete[] l_text;
	}

}


/// if kwdVal contains a hyphen it must end with "-CAR"
static void WarnCTYPE(fitsfile *fptr, const char* fileName, const char* kwdName)
{
char kwdVal[80];
int status = 0;
ReadFitsKey(fptr, kwdName, kwdVal, &status);
if (status==0) {
	const char* hyphen = strstr(kwdName, "-");
	if (hyphen && strcmp(hyphen, "-CAR"))
		cerr << "SkyMap Warning: In " << fileName << " " << kwdName << "=" << kwdName << endl;
	}
}

static void WarnCUNIT(fitsfile *fptr, const char* fileName, const char* kwdName)
{
char kwdVal[80];
int status = 0;
ReadFitsKey(fptr, kwdName, kwdVal, &status);
if (status==0 && strcmp(kwdVal, "deg"))
	cerr << "SkyMap Warning: In " << fileName << " " << kwdName << "=" << kwdVal << endl;
}


static void WarnCROTA(fitsfile *fptr, const char* fileName, const char* kwdName)
{
float kwdVal;
int status = 0;
ReadFitsKey(fptr, kwdName, &kwdVal, &status);
if (status==0 && kwdVal)
	cerr << "SkyMap Warning: In " << fileName << " " << kwdName << "=" << kwdVal << endl;
}



bool SkyMap::Load(const char* fileName)
{
fitsfile *fptr;
int status = 0;
fits_open_file(&fptr, fileName, READONLY, &status);
if (status)
	cerr << "SkyMap: failed opening file " << fileName << " Fits error: " << status << endl;
else {
	fits_movnam_hdu(fptr, IMAGE_HDU, const_cast<char*>("SKY"), 0, &status);
	if (status)
		cerr << "SkyMap: SKY HDU not found" << endl;
	else {
		int bitpix, naxis;
		long naxes[2];
		fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
		if (status)
			cerr << "SkyMap: Could not read the image" << endl;
		else if (naxis!=2)
			cerr << "SkyMap: Only 2D images are supported" << endl;
		if (status==0 && naxis==2) {
			int cols = naxes[0];
			int rows = naxes[1];
			float* values = new float[rows*cols];
			long fpixel[2] = { 1, 1 };
			fits_read_pix(fptr, TFLOAT, fpixel, rows*cols, 0, values, 0, &status);

			float scale, zero;
			float xRes, yRes, xRef, yRef, xPoint, yPoint;
			ReadFitsKey(fptr, "BSCALE", &scale, 1, &status);
			ReadFitsKey(fptr, "BZERO", &zero, 0, &status);
			ReadFitsKey(fptr, "CDELT1", &xRes, &status);
			ReadFitsKey(fptr, "CDELT2", &yRes, &status);
			ReadFitsKey(fptr, "CRPIX1", &xRef, &status);
			ReadFitsKey(fptr, "CRPIX2", &yRef, &status);
			ReadFitsKey(fptr, "CRVAL1", &xPoint, &status);
			ReadFitsKey(fptr, "CRVAL2", &yPoint, &status);

			if (status)
				cerr << "SkyMap: Error reading the basic keywords" << endl;
			else {
				WarnCTYPE(fptr, fileName, "CTYPE1");
				WarnCTYPE(fptr, fileName, "CTYPE2");
				WarnCUNIT(fptr, fileName, "CUNIT1");
				WarnCUNIT(fptr, fileName, "CUNIT2");
				WarnCROTA(fptr, fileName, "CROTA1");
				WarnCROTA(fptr, fileName, "CROTA2");
				}

			int kwrdCount;
			int* kwrdIndices;
			char* kwrdText;
			ReadUnknKeywords(fptr, kwrdCount, kwrdIndices, kwrdText, &status);
			/// PrintKeywords(fptr, &status);
			if (status)
				cerr << "SkyMap: Error reading ancillary keywords" << endl;
			else {
				delete[] m_values;
				m_values = values;
				m_rows = rows;
				m_cols = cols;
				SetRes(xRes, yRes);
				SetRefPoint(xRef, yRef, xPoint, yPoint);

				delete[] m_kwrdIndices;
				delete[] m_kwrdText;
				m_kwrdIndices = 0;
				m_kwrdText = 0;
				m_kwrdCount = kwrdCount;
				if (m_kwrdCount) {
					m_kwrdIndices = kwrdIndices;
					m_kwrdText = kwrdText;
					}
				}
			}
		}
	int closeStatus = 0;	/// Close this file in any case
	fits_close_file(fptr, &closeStatus);
	}
return status==0;
}


static void WriteStandardComments(fitsfile* fptr, int* status)
{
fits_write_comment(fptr, "FITS (Flexible Image Transport System) format is defined in 'Astronomy", status);
fits_write_comment(fptr, "and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H", status);
}


static void WriteUnknKeywords(fitsfile *fptr, int* status, int count, int* offsets, char* text)
{
for (int i=0; i<count; ++i) {
	char* card = text+offsets[i];
	fits_write_record(fptr, card, status);
	}
}



bool SkyMap::Save(const char* fileName) const
{
fitsfile* fptr;   /* FITS file pointer, defined in fitsio.h */
int status = 0;   /* CFITSIO status value MUST be initialized to zero! */

fits_create_file(&fptr, fileName, &status);
if (status) {
	cerr << "Create file failed" << endl;
	return false;
	}

time_t rawtime;
time(&rawtime);
struct tm* gmtTime = gmtime(&rawtime);
char isoTime[80];
strftime(isoTime, 80, "Written by PasteMap %Y-%m-%dT%H:%M:%S %Z", gmtTime);
char utcDate[16];
strftime(utcDate, 16, "%Y-%m-%d", gmtTime);

WriteKeyword(fptr, &status, "SIMPLE", true, isoTime);
WriteKeyword(fptr, &status, "BITPIX", 32, "Number of bits per data pixel");
WriteKeyword(fptr, &status, "NAXIS", 0, "Number of data axes");
WriteKeyword(fptr, &status, "EXTEND", true, "FITS data may contain extensions");
WriteKeyword(fptr, &status, "DATE", utcDate, "Creation UTC (CCCC-MM-DD) date of FITS header");
WriteStandardComments(fptr, &status);

long naxis = 2;
/// long naxes[2] = { m_rows, m_cols };
long naxes[2] = { m_cols, m_rows };
fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
long fpixel[2] = { 1, 1 };
fits_write_pix(fptr, TFLOAT, fpixel, m_rows*m_cols, m_values, &status);

	
/** Alrerady written by fits_write_pix():
WriteKeyword(fptr, &status, "XTENSION", "IMAGE", "IMAGE extension");
WriteKeyword(fptr, &status, "BITPIX", -32, "Number of bits per data pixel");
WriteKeyword(fptr, &status, "NAXIS", 2, "Number of data axes");
WriteKeyword(fptr, &status, "NAXIS1", m_rows, "Rows");
WriteKeyword(fptr, &status, "NAXIS2", m_cols, "Columns");
*/


WriteKeyword(fptr, &status, "CDELT1", m_xRes, "Horizontal resolution");
WriteKeyword(fptr, &status, "CRPIX1", m_xRef, "Reference column");
WriteKeyword(fptr, &status, "CRVAL1", m_xPoint, "Reference longitude");

WriteKeyword(fptr, &status, "CDELT2", m_yRes, "Vertical resolution");
WriteKeyword(fptr, &status, "CRPIX2", m_yRef, "Reference row");
WriteKeyword(fptr, &status, "CRVAL2", m_yPoint, "Reference latitude");

#ifdef MORE_KEYWORDS
WriteKeyword(fptr, &status, "CTYPE1", "GLON"); /// or GLON-CAR"
WriteKeyword(fptr, &status, "CUNIT1", "deg");
WriteKeyword(fptr, &status, "CTYPE2", "GLAT"); /// or "GLAT-CAR"
WriteKeyword(fptr, &status, "CUNIT2", "deg");

WriteKeyword(fptr, &status, "PIXCENT", true);
WriteKeyword(fptr, &status, "BUNIT", "INTENSITY");
WriteKeyword(fptr, &status, "PRIMTYPE", "RECT_MAP");
WriteKeyword(fptr, &status, "INSTRUME", "AGILE");
#endif

WriteKeyword(fptr, &status, "EXTNAME", "SKY");

WriteUnknKeywords(fptr, &status, m_kwrdCount, m_kwrdIndices, m_kwrdText);

int closeStatus = 0;	/// Close anyway
fits_close_file(fptr, &closeStatus);
return status==0;
}




SkyMap SkyMap::Split(int horSplit, int verSplit) const
{
if (horSplit==1 && verSplit==1)
	return *this;
SkyMap m(m_rows*verSplit, m_cols*horSplit, m_xRes/horSplit, m_yRes/verSplit, m_xPoint, m_yPoint, (m_xRef-1)*horSplit+1, (m_yRef-1)*verSplit+1);
m.CopyKeywords(*this);
m.Clear();
for (int i=0; i<m.m_rows; ++i)
	for (int j=0; j<m.m_cols; ++j)
		m.SetValue(i, j, Value(i/verSplit, j/horSplit));
return m;
}

SkyMap SkyMap::Collapse(int horBins, int verBins) const
{
if (horBins==1 && verBins==1)
	return *this;
SkyMap m(m_rows/verBins, m_cols/horBins, m_xRes*horBins, m_yRes*verBins, m_xPoint, m_yPoint, (m_xRef-1)/horBins+1, (m_yRef-1)/verBins+1);
m.CopyKeywords(*this);
m.Clear();
for (int i=0; i<m.m_rows; ++i)
	for (int j=0; j<m.m_cols; ++j) {
		float ave = 0;
		for (int ii=0; ii<verBins; ++ii)
			for (int jj=0; jj<horBins; ++jj)
				ave += Value(i*verBins+ii, j*horBins+jj);
		ave /= horBins*verBins;
		m.SetValue(i, j, ave);
		}
return m;
}




namespace
{

/// Zero ending array of prime numbers
	const int c_prime[] = {
		2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
		31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
		73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
		127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
		179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
		233, 239, 241, 251, 257, 263, 269, 271,   0 };


	/// a, b and c must be positive
	int CommonPrimes(int a, int b, int c)
	{
	int common = 1;
	int i = 0;
	while (c_prime[i] && c_prime[i]<=a && c_prime[i]<=b) {
		int p = c_prime[i];
		if (a%p==0 && b%p==0 && c%p==0) {
			common *= p;
			a /= p;
			b /= p;
			c /= p;
			}
		else
			++i;
		}
	return common;
	}

}


extern SkyMap PasteMap(const SkyMap& map1, const SkyMap& map2, float binSize)
{
int xbin1 = abs(int(map1.XRes()*100));
int ybin1 = abs(int(map1.YRes()*100));
int xbin2 = abs(int(map2.XRes()*100));
int ybin2 = abs(int(map2.YRes()*100));
int xbin3 = binSize<=0 ? xbin1 : int(binSize*100);
int ybin3 = binSize<=0 ? ybin1 : int(binSize*100);

int xbin = CommonPrimes(xbin1, xbin2, xbin3);
int ybin = CommonPrimes(ybin1, ybin2, ybin3);

SkyMap map = map1.Split(ybin1/ybin, xbin1/xbin);

float x, y;
int row, col;
for (int i=0; i<map.Rows(); ++i)
	for (int j=0; j<map.Cols(); ++j) {
		map.Bin2Point(i, j, &x, &y);
		if (map2.Point2Bin(x, y, &row, &col))
			map(i, j) = map2(row, col);
		}
if (binSize>0)
	map = map.Collapse(xbin3/xbin, ybin3/ybin);
return map;
}
