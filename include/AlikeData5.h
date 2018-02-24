/*
 * Copyright (c) 2011
 *     Tomaso Contessi (Nuove Idee sas)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 * GPL License
*/

#ifndef _ALIKE_DATA5_
#define _ALIKE_DATA5_


#include <iostream>
#include <vector>
#include <string>
#include "AgileMap.h"

using std::vector;
using std::string;


/// Loading a file with this name will return an empty data structure
#define c_missingFileName "None"


/// FixFlag values
typedef int FixFlag;
const FixFlag AllFixed  = 0;
const FixFlag FluxFree  = 1;
const FixFlag PosFree   = 2;
const FixFlag IndexFree = 4;
const FixFlag Par2Free =  8;
const FixFlag Par3Free = 16;
/// typedef enum { AllFixed=0, FluxFree = 1, PosFree = 2, IndexFree = 4 } FixFlag;


/// //////////////////////////////////////////////////////////////////
///
/// Class for the exposure correction
///
/*
class ExpCorr
{
public:
	ExpCorr(): m_thetaArr(0), m_corrArr(0), m_count(0) {}
	ExpCorr(const char* fileName): m_thetaArr(0), m_corrArr(0), m_count(0) { Load(fileName); }
	ExpCorr(const ExpCorr& expCorr) { Copy(expCorr); }
	~ExpCorr() { Destroy(); }
	bool Load(const char* fileName);

	ExpCorr& operator=(const ExpCorr& expCorr);

	float GetCorrection(float theta) const;
private:	/// Data
	float* m_thetaArr;
	float* m_corrArr;
	long   m_count;
	void Destroy() { delete[] m_thetaArr; delete[] m_corrArr; }
	void Clear() { m_count=0; m_thetaArr=m_corrArr=0; }
	void Copy(const ExpCorr& expCorr);
};
*/


/// File format: text file each line containing the following space separated values
/// flux srcL srcB index fixflag minTS label [loclimit]
/// A line beginning with an exclamation mark is a comment
/// The label is a string contain no spaces, all the other items are number
/// An absent or zero loclimit means the source location is not limited
struct SourceData
{
/// Input data
	string  label;
	FixFlag fixflag;
	double  minTS;
	double  loclimit;	/// sources can move in [l+loclimit, l-loclimit]

/// Input/output data
	double flux;
	double srcL;
	double srcB;
	double index;
	//0 - PL (k E^-{\index})
	//1 - PLExpCutoff k E^-{\index} e^ ( - E / E_c ) -> par2 = E_c
	int typefun; //0 PowerLaw - 1 PLExpCutoff - 2 LogParabola
	double par2;
	double par3;

/// Output data
	double TS;
	double gal;
	double iso;
	double fluxul;

	SourceData(): label(), fixflag(0), minTS(0), loclimit(0), flux(0), srcL(0), srcB(0), index(0), typefun(0), par2(0), par3(0), TS(0), gal(0), iso(0), fluxul(0) {}
	SourceData(const SourceData& another):
		label(another.label), fixflag(another.fixflag), minTS(another.minTS), loclimit(another.loclimit), flux(another.flux), srcL(another.srcL), srcB(another.srcB), index(another.index), typefun(another.typefun), par2(another.par2), par3(another.par3), TS(another.TS), gal(another.gal), iso(another.iso), fluxul(another.fluxul) {}
	
	/// void GetResultsFrom(const AlikeSourceMap& map);
	
	void Print(std::ostream& oFile) const;
	void VerbosePrint(std::ostream& oFile) const;
};



class SourceDataArray
{
public:
	SourceDataArray(int length=0): m_dataArr(length) {}
	SourceDataArray(const SourceDataArray& another): m_dataArr(another.m_dataArr) {}

	friend SourceDataArray operator+(const SourceDataArray& arr1, const SourceDataArray& arr2);

	int Count() const { return m_dataArr.size(); }
	void Append(const SourceData& data) { m_dataArr.push_back(data); }
	void SortByFlux();
	const SourceData& operator[](int index) const { return m_dataArr[index]; }
	SourceData& operator[](int index) { return m_dataArr[index]; }
	const SourceData& operator[](string label) const;
	SourceData& operator[](string label);
	
	void Print(std::ostream& oFile) const;
private:
	vector<SourceData> m_dataArr;
	int IndexOf(string label) const;
};


extern SourceDataArray ReadSourceFile(const char* srcfilename);







/// File format: text file each line containing the following space separated values
/// CtsName ExpName GasName Theta GalCoeff IsoCoeff
/// A line beginning with an exclamation mark is a comment
/// A negative value of GalCoeff or IsoCoeff means a variable coefficient

class MapList
{
public:
	MapList(): m_count(0) {}
	~MapList() {}
	bool Read(const char* listFileName);

	/// Maps
	int Count() const { return m_count; }
	const char* CtsName(int i) const { return m_ctslist[i].c_str(); }
	const char* ExpName(int i) const { return m_explist[i].c_str(); }
	const char* GasName(int i) const { return m_gaslist[i].c_str(); }
	double Theta(int i) const { return m_thetaArr[i]; }
	double GalCoeff(int i) const { return m_galCoeffArr[i]; }
	double IsoCoeff(int i) const { return m_isoCoeffArr[i]; }
private:
	int m_count;
	vector<string> m_ctslist;
	vector<string> m_explist;
	vector<string> m_gaslist;
	vector<double> m_thetaArr;
	vector<double> m_galCoeffArr;
	vector<double> m_isoCoeffArr;
};


class MapCoeff
{
public:
	MapCoeff() { Clear(); }
	MapCoeff(double theta, double galCoeff, double isoCoeff);
	MapCoeff(const MapCoeff& another) { Copy(another); }
	~MapCoeff() { Destroy(); }
	MapCoeff& operator=(const MapCoeff& another) { if (this!=&another) { Destroy(); Copy(another); } return *this; }
	bool Load(const MapList& maplist);

	int Length() const { return m_length; }

	double Theta(int i) const { return m_thetaArr[i]; }
	double GalCoeff(int i) const { return m_galCoeffs[i]; }
	double IsoCoeff(int i) const { return m_isoCoeffs[i]; }
protected:
	int       m_length;
	double*   m_thetaArr;
	double*   m_galCoeffs;
	double*   m_isoCoeffs;
	void Destroy();
	void Clear();
	void Copy(const MapCoeff& another);
};




class MapMaps
{
public:
	MapMaps() { Clear(); }
	MapMaps(const AgileMap& ctsMap, const AgileMap& expMap, const AgileMap& gasMap);
	MapMaps(const MapMaps& another) { Copy(another); }
	~MapMaps() { Destroy(); }
	MapMaps& operator=(const MapMaps& another) { if (this!=&another) { Destroy(); Copy(another); } return *this; }
	bool Load(const MapList& maplist, bool skipCts=false);

	bool SetCtsMap(const AgileMap& ctsMap, int i);	/// Deprecated
	bool ReplaceCtsMaps(AgileMap* ctsMapArr);			/// Array of Count() maps will be owned by this class

	int Count() const { return m_mapCount; }

	const AgileMap& CtsMap(int i) const { return m_ctsMaps[i]; }
	const AgileMap& ExpMap(int i) const { return m_expMaps[i]; }
	const AgileMap& GasMap(int i) const { return m_gasMaps[i]; }

	double MapCounts(int i) const { return m_countsArr[i]; }

	bool Aligned() const { return m_alignedMaps; }
	double EnergyInf() const { return m_energyInf; }
	double EnergySup() const { return m_energySup; }

protected:
	int       m_mapCount;
	AgileMap* m_ctsMaps;
	AgileMap* m_expMaps;
	AgileMap* m_gasMaps;
	int*      m_loadFlags;	/// 0=OK; 1=No cts; 2=ctsErr; 4=expErr; 8=gasErr; 16 = cts-exp mismatch; 32 = exp-gas mismatch
	int*      m_countsArr;
	double    m_energyInf;
	double    m_energySup;
	bool      m_alignedMaps;
	void Destroy();
	void Clear();
	void Copy(const MapMaps& another);
};


class MapData: public MapMaps, public MapCoeff
{
public:
	MapData(): MapMaps(), MapCoeff() {}
	MapData(const AgileMap& ctsMap, const AgileMap& expMap, const AgileMap& gasMap, double theta, double galCoeff, double isoCoeff): MapMaps(ctsMap, expMap, gasMap), MapCoeff(theta, galCoeff, isoCoeff) {}
	MapData(const MapData& another): MapMaps(another), MapCoeff(another) {}
	~MapData() {}
	MapData& operator=(const MapData& another) { if (this!=&another) { Destroy(); Copy(another); } return *this; }
	bool Load(const MapList& maplist, bool skipCts=false);
	void Print() const;

private:
	void Destroy() { MapMaps::Destroy(); MapCoeff::Destroy(); }
	void Clear() { MapMaps::Clear(); MapCoeff::Clear(); }
	void Copy(const MapData& another) { MapMaps::Copy(another); MapCoeff::Copy(another); }
};




/// ////////////////////////////////////////////////////////////////////////
///
/// Extended Sources
///



class ExtList
{
public:
	ExtList(): m_mapCount(0), m_extCount(0) {}
	~ExtList() {}
	bool Read(const char* listFileName);

	int MapCount() const { return m_mapCount; }
	int ExtCount() const { return m_extCount; }

	const char* ExtName(int i) const { return m_extNames[i].c_str(); } 
	double ExtFlux(int i) const { return m_fluxArr[i]; } 
	double ExtIndex(int i) const { return m_indexArr[i]; } 
	const char* ExtFileName(int map, int i) const { return m_extMaps[Idx(map, i)].c_str(); } 
private:
	int m_mapCount;
	int m_extCount;
	vector<string> m_extNames;
	vector<double> m_fluxArr;
	vector<double> m_indexArr;
	vector<string> m_extMaps;
	void Clear() { m_mapCount=m_extCount=0; m_extNames.clear(); m_fluxArr.clear(); m_indexArr.clear(); m_extMaps.clear(); }
	/// int Idx(int map, int i) const { return i*m_mapCount+map; }
	int Idx(int map, int i) const { return map*m_extCount+i; }
};


class ExtData
{
public:
	ExtData() { Clear(); }
	ExtData(const ExtData& another) { Copy(another); }
	~ExtData() { Destroy(); }

	ExtData& operator=(const ExtData& another) { if (this!=&another) { Destroy(); Copy(another); } return *this; }
	ExtData& operator=(const ExtList& extlist) { Destroy(); Load(extlist); return *this; }

	bool Load(const ExtList& extlist);

	int MapCount() const { return m_mapCount; }
	int ExtCount() const { return m_extCount; }

	const AgileMap& ExtMap(int map, int i) const { return m_extArr[Idx(map, i)]; }
	AgileMap& ExtMap(int map, int i) { return m_extArr[Idx(map, i)]; }

	double Flux(int i) const { return m_fluxArr[i]; }
	double Index(int i) const { return m_indexArr[i]; }
	const char* Name(int i) const { return &m_nameBuff[m_nameOfs[i]]; }

private:
	int       m_mapCount;
	int       m_extCount;
	AgileMap* m_extArr;
	double*   m_fluxArr;
	double*   m_indexArr;
	char*     m_nameBuff;
	int*      m_nameOfs;
	int       m_nameBuffSize;
	vector<string> m_extNames;
	/// int Idx(int map, int i) const { return i*m_mapCount+map; }
	int Idx(int map, int i) const { return map*m_extCount+i; }
	void Destroy() { delete[] m_extArr; delete[] m_fluxArr; delete[] m_indexArr; delete[] m_nameBuff; delete[] m_nameOfs; }
	void Clear() { m_mapCount=m_extCount=m_nameBuffSize=0; m_extArr=0; m_fluxArr=m_indexArr=0; m_nameBuff=0; m_nameOfs=0; }
	void Copy(const ExtData& another);
};

#endif
