////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       2011-03-28
//       Authors: Tomaso Contessi, IASF-Milano
//
// INFO
//     Class for handlig an image representing a rectangular SKYMAP.
//     Methods for loading and saving to a FITS file
//     The FITS file shall contain the keywords:
//     CDELT1, CDELT2, CRPIX1, CRPIX2, CRVAL1, CRVAL2, and optionally BSCALE and BZERO
//     All the other keywords possibly present in the FITS file are loaded, 
//     preserved in all the copies, and saved.
//
// HISTORY
//      Version 1.0
//      Date 2011-03-28
//      Author Tomaso Contessi
//
//      Version 1.1
//      Date 2013-07-27
//      Addressed an issue related to BSCALE and BZERO
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




#ifndef _SKY_MAP_
#define _SKY_MAP_



const char* const stdKeywords[] = {
"CTYPE1  = 'GLON-CAR'           /",
"CUNIT1  = 'deg     '           /",
"CTYPE2  = 'GLAT-CAR'           /",
"CUNIT2  = 'deg     '           /",
"PIXCENT =                    T /",
"BUNIT   = 'INTENSITY'          /",
"PRIMTYPE= 'RECT_MAP'           /",
"INSTRUME= 'AGILE   '           /",
0 };

const char* const stdKeywords_noCAR[] = {
"CTYPE1  = 'GLON    '           /",
"CUNIT1  = 'deg     '           /",
"CTYPE2  = 'GLAT    '           /",
"CUNIT2  = 'deg     '           /",
"PIXCENT =                    T /",
"BUNIT   = 'INTENSITY'          /",
"PRIMTYPE= 'RECT_MAP'           /",
"INSTRUME= 'AGILE   '           /",
0 };




class SkyMap
{
public:
	SkyMap() { Reset(); SetDefaults(); }
	SkyMap(int rows, int cols) { Reset(); Resize(rows, cols); SetDefaults(); }
	SkyMap(int rows, int cols, float xRes, float yRes, float xPoint, float yPoint, float xRef=-1, float yRef=-1);
	SkyMap(const SkyMap& other) { Reset(); Copy(other); }
	~SkyMap() { Destroy(); }

	SkyMap& operator=(const SkyMap& other) { if (this!=&other) Copy(other); return *this; }

	void SetKeywords(const char* const keywords[]);

	bool Load(const char* fileName);
	bool Save(const char* fileName) const;
	void Print() const;

	SkyMap Split(int horFactor, int verFactor) const; /// Split each bin into a rect of bins (same value)
	SkyMap Collapse(int horBins, int verBins) const;  /// Collapse a rect into one bin (mean value). Rows and cols get cut

/*
	void Clear(float scale=1, float zero=0);   /// Clear all the physical values
	void Scale(float newScale, float newZero); /// Change scale and zero keeping the physical values
*/
	int  Rows() const { return m_rows; }
	int  Cols() const { return m_cols; }
	void Range(float* minX, float* maxX, float* minY, float* maxY) const; /// Any param can be zero


	float& operator()(int row, int col) { return m_values[row*m_cols+col]; }
	float operator()(int row, int col) const { return m_values[row*m_cols+col]; }

	void Clear() { SetAllValues(0.0f); }
	void SetAllValues(float value);

	float& Value(int row, int col) { return m_values[row*m_cols+col]; }
	float Value(int row, int col) const { return m_values[row*m_cols+col]; }
	void SetValue(int row, int col, float value) { Value(row, col) = value; }

	bool Point2Bin(float x, float y, int* row, int* col) const; /// Convert coords to bin. Any param can be zero
	bool Bin2Point(int row, int col, float* x, float* y) const; /// Convert bin to coords. Any param can be zero

	float XRes() const { return m_xRes; }
	float YRes() const { return m_yRes; }
	void  SetRes(float xRes, float yRes) { m_xRes = xRes; m_yRes = yRes; }

	void GetRefPoint(float* xRef, float* yRef, float* xPoint, float* yPoint) const; /// Any param can be zero
	void SetRefPoint(float xRef, float yRef, float xPoint, float yPoint);

private:
	float* m_values;
	int    m_rows;
	int    m_cols;

	float  m_xRes;		/// CDELT1 Horizontal degrees per bin
	float  m_yRes;		/// CDELT2 Vertical degrees per bin
	float  m_xRef;		/// CRPIX1 column of the reference bin
	float  m_yRef;		/// CRPIX2 row of the reference bin
	float  m_xPoint;	/// CRVAL1 x coord pointed by m_xRef
	float  m_yPoint;	/// CRVAL2 y coord pointed by m_yRef

	/// Ancillary data
	float m_minX;
	float m_maxX;
	float m_minY;
	float m_maxY;
	void SetLimits();

	void SetDefaults() { m_xRes=m_yRes=1; m_xPoint=m_yPoint=m_xRef=m_yRef=m_minX=m_maxX=m_minY=m_maxY=0; }

	/// Unknown keywords
	int   m_kwrdCount;
	int*  m_kwrdIndices;
	char* m_kwrdText;
	void CopyKeywords(const SkyMap& other);

	/// Internal operations
	void Reset() { m_values=0; m_rows=m_cols=0; m_kwrdCount=0; m_kwrdIndices=0; m_kwrdText=0; }
	void Destroy() { delete[] m_values; delete[] m_kwrdIndices; delete[] m_kwrdText; }
	void Resize(int rows, int cols);
	void Copy(const SkyMap& other);
};



/**
Paste map1 on map2, giving map3.
The width of the maps bin must be multiple of 0.01 degrees per bin, and the resulting map will have the resolution given by binSize, or the same of map1 if binSize <=0.
*/


SkyMap PasteMap(const SkyMap& map1, const SkyMap& map3, float binSize=-1);



#endif
