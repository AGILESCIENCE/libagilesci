



#ifndef _FITS_UTILS_
#define _FITS_UTILS_


#include "fitsio.h"


/// Functions reading and writing keys (they do nothing if status already points to an error)
/// They return true if the keyword was actually read or written

/// Reading optional keyword. They don't set the error KEY_NO_EXIST, and assign the default vsalue in that case
bool ReadFitsKey(fitsfile* fptr, const char* name, float* value, float defaultVal, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, double* value, double defaultVal, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, int* value, int defaultVal, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, long* value, long defaultVal, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, char* value, const char* defaultVal, int* status);

/// Reading mandatory keyword.
bool ReadFitsKey(fitsfile* fptr, const char* name, float* value, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, double* value, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, int* value, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, long* value, int* status);
bool ReadFitsKey(fitsfile* fptr, const char* name, char* value, int* status);

/// Writing keywords. They write it again if it was already present
bool WriteFitsKey(fitsfile *fptr, const char* name, bool value, int* status, const char* comment=0);
bool WriteFitsKey(fitsfile *fptr, const char* name, float value, int* status, const char* comment=0);
bool WriteFitsKey(fitsfile *fptr, const char* name, double value, int* status, const char* comment=0);
bool WriteFitsKey(fitsfile *fptr, const char* name, int value, int* status, const char* comment=0);
bool WriteFitsKey(fitsfile *fptr, const char* name, long value, int* status, const char* comment=0);
bool WriteFitsKey(fitsfile *fptr, const char* name, const char* value, int* status, const char* comment=0);


bool WriteFitsAstroComments(fitsfile* fptr, int* status);
bool WriteFitsDate(fitsfile* fptr, int* status); /// Write the DATE keyword with UTC time


class FitsFile
{
public:
	FitsFile(): m_status(0), m_fptr(0) {}
	FitsFile(fitsfile* fptr, int status=0): m_status(status), m_fptr(fptr) {}
	FitsFile(const char* fileName, int iomode=READONLY): m_status(0), m_fptr(0) { Open(fileName, iomode); }
	~FitsFile() { Close(); }

	operator fitsfile*() { return m_fptr; }
	operator int*() { return &m_status; }

	fitsfile* File() { return m_fptr; }
	int Status() const { return m_status; }
	bool Ok() const { return !m_status; }
	void Error(char* text) const { fits_get_errstatus(m_status, text); } // Pass a 32 char string

	void Detach() { m_fptr=0; m_status=0; }
	void Attach(fitsfile* fptr, int status=0) { m_fptr=fptr; m_status=status; } /// Warning, it detaches the previous file
	void SetStatus(int status) { m_status = status; }

	bool MoveAbsHDU(int hdu)
		{ if (!m_fptr) return false; fits_movabs_hdu(m_fptr, hdu, NULL, &m_status); return !m_status; }
	bool CreateHDU()
		{ if (!m_fptr) return false; fits_create_hdu(m_fptr, &m_status); return !m_status; }

	bool Open(const char* fileName, int iomode=READONLY) /// READWRITE otherwise
		{ if (m_fptr) return false; fits_open_file(&m_fptr, const_cast<char*>(fileName), iomode, &m_status); return !m_status; }
	bool Create(const char* fileName)
		{ if (m_fptr) return false; fits_create_file(&m_fptr, const_cast<char*>(fileName), &m_status); return !m_status; }

	/// Close and Delete functions keep the error of the last operation performed before closing or deleting, but they try to close or delete anyway
	bool Close()
		{ if (!m_fptr) return false; int status=0; fits_close_file(m_fptr, &status); m_fptr=0; return !status; }
	bool Delete()
		{ if (!m_fptr) return false; int status=0; fits_delete_file(m_fptr, &status); m_fptr=0; return !status; }


	bool ReadKey(const char* name, float* value, float defaultVal)
		{ return ::ReadFitsKey(m_fptr, name, value, defaultVal, &m_status); }
	bool ReadKey(const char* name, double* value, double defaultVal)
		{ return ::ReadFitsKey(m_fptr, name, value, defaultVal, &m_status); }
	bool ReadKey(const char* name, long* value, long defaultVal)
		{ return ::ReadFitsKey(m_fptr, name, value, defaultVal, &m_status); }
	bool ReadKey(const char* name, char* value, const char* defaultVal)
		{ return ::ReadFitsKey(m_fptr, name, value, defaultVal, &m_status); }

	bool ReadKey(const char* name, float* value) { return ::ReadFitsKey(m_fptr, name, value, &m_status); }
	bool ReadKey(const char* name, double* value) { return ::ReadFitsKey(m_fptr, name, value, &m_status); }
	bool ReadKey(const char* name, long* value) { return ::ReadFitsKey(m_fptr, name, value, &m_status); }
	bool ReadKey(const char* name, char* value) { return ::ReadFitsKey(m_fptr, name, value, &m_status); }

	bool WriteKey(const char* name, bool value, const char* comment=0)
		{ return ::WriteFitsKey(m_fptr, name, value, &m_status, comment); }
	bool WriteKey(const char* name, float value, const char* comment=0)
		{ return ::WriteFitsKey(m_fptr, name, value, &m_status, comment); }
	bool WriteKey(const char* name, double value, const char* comment=0)
		{ return ::WriteFitsKey(m_fptr, name, value, &m_status, comment); }
	bool WriteKey(const char* name, int value, const char* comment=0)
		{ return ::WriteFitsKey(m_fptr, name, value, &m_status, comment); }
	bool WriteKey(const char* name, long value, const char* comment=0)
		{ return ::WriteFitsKey(m_fptr, name, value, &m_status, comment); }
	bool WriteKey(const char* name, const char* value, const char* comment=0)
		{ return ::WriteFitsKey(m_fptr, name, value, &m_status, comment); }

	bool WriteAstroComments() { return WriteFitsAstroComments(m_fptr, &m_status); }
	bool WriteDate() { return WriteFitsDate(m_fptr, &m_status); }

	bool ColumnInfo(const char* colName, int* colnum, long* nrows);

	/// bool ReadArray(const char* colName, VecF& floatArr);
	/// bool ReadArray(const char* colName, VecD& doubleArr);

	/// bool GetImageSize(MatrixOf<long>& naxes);

private:
	int       m_status;
	fitsfile* m_fptr;

private: /// No implementation
	FitsFile(const FitsFile&);
	const FitsFile& operator=(const FitsFile&);
};


#endif
