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
#include <fstream>
#include "FitsUtils.h"

void CopyFile(const char* iname, const char* oname)
{
    std::ifstream ifs(iname, std::ios::binary);
    std::ofstream ofs(oname, std::ios::binary);
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
}

bool ReadFitsKey(fitsfile* fptr, const char* name, float* value, int* status)
{
if (!*status)
	fits_read_key(fptr, TFLOAT, const_cast<char*>(name), value, NULL, status);
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, double* value, int* status)
{
if (!*status)
	fits_read_key(fptr, TDOUBLE, const_cast<char*>(name), value, NULL, status);
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, int* value, int* status)
{
if (!*status)
	fits_read_key(fptr, TINT, const_cast<char*>(name), value, NULL, status);
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, long* value, int* status)
{
if (!*status)
	fits_read_key(fptr, TLONG, const_cast<char*>(name), value, NULL, status);
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, char* value, int* status)
{
if (!*status)
	fits_read_key(fptr, TSTRING, const_cast<char*>(name), value, NULL, status);
return !*status;
}



bool ReadFitsKey(fitsfile* fptr, const char* name, float* value, float defaultVal, int* status)
{
if (!*status) {
	ReadFitsKey(fptr, name, value, status);
	if (KEY_NO_EXIST==*status) {
		*value = defaultVal;
		*status = 0;
		return false;
		}
	}
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, double* value, double defaultVal, int* status)
{
if (!*status) {
	ReadFitsKey(fptr, name, value, status);
	if (KEY_NO_EXIST==*status) {
		*value = defaultVal;
		*status = 0;
		return false;
		}
	}
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, int* value, int defaultVal, int* status)
{
if (!*status) {
	ReadFitsKey(fptr, name, value, status);
	if (KEY_NO_EXIST==*status) {
		*value = defaultVal;
		*status = 0;
		return false;
		}
	}
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, long* value, long defaultVal, int* status)
{
if (!*status) {
	ReadFitsKey(fptr, name, value, status);
	if (KEY_NO_EXIST==*status) {
		*value = defaultVal;
		*status = 0;
		return false;
		}
	}
return !*status;
}

bool ReadFitsKey(fitsfile* fptr, const char* name, char* value, const char* defaultVal, int* status)
{
if (!*status) {
	ReadFitsKey(fptr, name, value, status);
	if (KEY_NO_EXIST==*status) {
		strcpy(value, defaultVal);
		*status = 0;
		return false;
		}
	}
return !*status;
}


bool WriteFitsKey(fitsfile *fptr, const char* name, bool value, int* status, const char* comment)
{
if (!*status)
	fits_write_key(fptr, TLOGICAL, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool WriteFitsKey(fitsfile *fptr, const char* name, float value, int* status, const char* comment)
{
if (!*status)
	fits_write_key(fptr, TFLOAT, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool WriteFitsKey(fitsfile *fptr, const char* name, double value, int* status, const char* comment)
{
if (!*status)
	fits_write_key(fptr, TDOUBLE, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool WriteFitsKey(fitsfile *fptr, const char* name, int value, int* status, const char* comment)
{
if (!*status)
	fits_write_key(fptr, TINT, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool WriteFitsKey(fitsfile *fptr, const char* name, long value, int* status, const char* comment)
{
if (!*status)
	fits_write_key(fptr, TLONG, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool WriteFitsKey(fitsfile *fptr, const char* name, const char* value, int* status, const char* comment)
{
if (!*status)
	fits_write_key(fptr, TSTRING, const_cast<char*>(name), const_cast<char*>(value), const_cast<char*>(comment), status);
return *status;
}


bool UpdateFitsKey(fitsfile *fptr, const char* name, bool value, int* status, const char* comment)
{
if (!*status)
	fits_update_key(fptr, TLOGICAL, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool UpdateFitsKey(fitsfile *fptr, const char* name, float value, int* status, const char* comment)
{
if (!*status)
	fits_update_key(fptr, TFLOAT, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool UpdateFitsKey(fitsfile *fptr, const char* name, double value, int* status, const char* comment)
{
if (!*status)
	fits_update_key(fptr, TDOUBLE, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool UpdateFitsKey(fitsfile *fptr, const char* name, int value, int* status, const char* comment)
{
if (!*status)
	fits_update_key(fptr, TINT, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool UpdateFitsKey(fitsfile *fptr, const char* name, long value, int* status, const char* comment)
{
if (!*status)
	fits_update_key(fptr, TLONG, const_cast<char*>(name), &value, const_cast<char*>(comment), status);
return *status;
}

bool UpdateFitsKey(fitsfile *fptr, const char* name, const char* value, int* status, const char* comment)
{
if (!*status)
	fits_update_key(fptr, TSTRING, const_cast<char*>(name), const_cast<char*>(value), const_cast<char*>(comment), status);
return *status;
}



bool FitsFile::ColumnInfo(const char* colName, int* colnum, long* nrows)
{
if (m_status)
	return false;
fits_get_colnum(m_fptr, CASEINSEN, const_cast<char*>(colName), colnum, &m_status);
int typecode = 0;
long width = 0;
fits_get_eqcoltype(m_fptr, *colnum, &typecode, nrows, &width, &m_status);
return !m_status;
}




/** If we want to include Matrix.h

bool FitsFile::GetImageSize(MatrixOf<long>& naxes)
{
int naxis = 0;
fits_get_img_dim(m_fptr, &naxis, &m_status);
if (!m_status) {
	naxes.ReshapeTo(1, naxis);
	fits_get_img_size(m_fptr, naxis, naxes, &m_status);
	}
return !m_status;
}



bool FitsFile::ReadArray(const char* colName, VecF& floatArr)
{
long nrows = 0;
int  colnum = 0;
if (ArrayLength(colName, &colnum, &nrows) {
	floatArr.ReshapeTo(1, nrows);
	fits_read_col(m_fptr, TFLOAT, colnum, 1, 1, nrows, NULL, floatArr, NULL, &m_status);
	}
return !m_status;
}

bool FitsFile::ReadArray(const char* colName, VecD& doubleArr)
{
long nrows = 0;
int  colnum = 0;
if (ArrayLength(colName, &colnum, &nrows) {
	floatArr.ReshapeTo(1, nrows);
	fits_read_col(m_fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, floatArr, NULL, &m_status);
	}
return !m_status;
}

*/


/**
extern bool ReadOptFitsKey(fitsfile* fptr, const char* name, double* value, double defaultVal, int* status)
{
if (ReadFitsKey(fptr, name, value, status))
	return true;
if (KEY_NO_EXIST==*status) {
	*status = 0;
	
*value = defaultVal;
if (!*status) {
	int status2 = 0;
	double value2;
	fits_read_key(fptr, TDOUBLE, const_cast<char*>(name), &value2, NULL, &status2);
	if (status2==0) {
		*value = value2;
		return true;
		}
	else if (status2!=KEY_NO_EXIST)
		*status = status2;
	}
return false;
}

extern bool ReadOptFitsKey(fitsfile* fptr, const char* name, long* value, long defaultVal, int* status)
{
*value = defaultVal;
if (!*status) {
	int status2 = 0;
	double value2;
	fits_read_key(fptr, TLONG, const_cast<char*>(name), &value2, NULL, &status2);
	if (status2==0) {
		*value = value2;
		return true;
		}
	else if (status2!=KEY_NO_EXIST)
		*status = status2;
	}
return false;
}

extern bool ReadOptFitsKey(fitsfile* fptr, const char* name, char* value, const char* defaultVal, int* status)
{
strcpy(value, defaultVal);
if (!*status) {
	int status2 = 0;
	char value2[1024];
	fits_read_key(fptr, TSTRING, const_cast<char*>(name), value2, NULL, &status2);
	if (status2==0) {
		strcpy(value, value2);
		return true;
		}
	else if (status2!=KEY_NO_EXIST)
		*status = status2;
	}
return false;
}
*/
