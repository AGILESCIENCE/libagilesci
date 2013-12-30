


#include <cstring>
#include "FitsUtils.h"



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