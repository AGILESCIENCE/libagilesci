



#include <cstring>

#include <sstream>
#include <fstream>

#include <algorithm>

#include "FitsUtils.h"
#include "AlikeData5.h"


using namespace std;


/// //////////////////////////////////////////////////////////////////
///
/// Class for the exposure correction
///

/*

void ExpCorr::Copy(const ExpCorr& expCorr)
{
m_count = expCorr.m_count;
if (!m_count)
	Clear();
else {
	m_thetaArr = new float[m_count];
	m_corrArr = new float[m_count];
	for (int i=0; i<m_count; ++i) {
		m_thetaArr[i] = expCorr.m_thetaArr[i];
		m_corrArr[i] = expCorr.m_corrArr[i];
		}
	}
}


ExpCorr& ExpCorr::operator=(const ExpCorr& expCorr)
{
if (this!=&expCorr) {
	Destroy();
	Copy(expCorr);
	}
return *this;
}


bool ExpCorr::Load(const char* fileName)
{
Destroy();
Clear();
if (fileName==0)
	return true;
if (fileName[0]==0)
	return true;
if (strcmp(fileName, c_missingFileName)==0)
	return true;
///	if (fileName==0 || fileName[0]==0 || strcmp(fileName, c_missingFileName)==0)
/// 		return true;
int status = 0;
fitsfile* fptr;
status = fits_open_table(&fptr, fileName, READONLY, &status);
if (status) {
	cout << "Information: file " << fileName << " not found" << endl;
	return false;
	}
status = fits_get_num_rows(fptr, &m_count, &status);
if (status) {
	cout << "ERROR accessing file " << fileName << endl;
	fits_close_file(fptr, &status);
	return false;
	}
m_thetaArr = new float[m_count];
m_corrArr = new float[m_count];
status = fits_read_col(fptr, TFLOAT, 1, 1, 1, m_count, 0, m_thetaArr, 0, &status);
status = fits_read_col(fptr, TFLOAT, 2, 1, 1, m_count, 0, m_corrArr, 0, &status);
if (status) {
	cout << "ERROR reading file " << fileName << endl;
	delete[] m_thetaArr;
	delete[] m_corrArr;
	m_thetaArr = 0;
	m_corrArr = 0;
	m_count = 0;
	}
// cout << m_count << ": " ;
// for (int i=0; i<m_count; i++) {cout << m_thetaArr[i] << "," << m_corrArr[i] << "; " ;}
// cout << endl;
fits_close_file(fptr, &status);
return status==0;
}

float ExpCorr::GetCorrection(float theta) const
{
if (!m_count)
	return 1;
if (theta<=m_thetaArr[0] || m_count==1)
	return m_corrArr[0];
if (theta>=m_thetaArr[m_count-1])
	return m_corrArr[m_count-1];
int i = 0;
for (; theta>m_thetaArr[i]; ++i) ;
 --i;
float deltaTheta = m_thetaArr[i+1]-m_thetaArr[i];
if (deltaTheta==0)
	return m_corrArr[i];
float deltaCorr = m_corrArr[i+1]-m_corrArr[i];
float fraction = (theta-m_thetaArr[i])/deltaTheta;
// cout << theta << "," << m_corrArr[i]+deltaCorr*fraction << endl;
// cout << i << " " << m_thetaArr[i] << " " << m_corrArr[i] << " " << deltaCorr << endl;
return m_corrArr[i]+deltaCorr*fraction;
}
*/


/// /////////////////////////////////////////////////////////
///
/// SourceData and SourceDataArray
///
/// Construction parameters for AlikeSourceMap
///


void SourceData::Print(ostream& oFile) const
{
oFile <<flux << " " << srcL << " " << srcB << " " << index << " " << fixflag << " " << sqrt(minTS) << " " << label;
if (loclimit!=0.0)
	oFile << " " << loclimit;
oFile << endl;
}

void SourceData::VerbosePrint(ostream& oFile) const
{
oFile << label << ": Flux = " << flux << ", L = " << srcL << ", B = " << srcB << ", Index = " << index << ", FixFlag = " << fixflag << ", sqrt(minTS) = " << sqrt(minTS);
if (loclimit!=0.0)
	oFile << ", loclimit = " << loclimit;
oFile << endl;
}


SourceDataArray operator+(const SourceDataArray& arr1, const SourceDataArray& arr2)
{
int l1 = arr1.Count();
int l2 = arr2.Count();
SourceDataArray arr(l1+l2);
for (int i=0; i<l1; ++i)
	arr[i] = arr1[i];
for (int i=0; i<l2; ++i)
	arr[i+l1] = arr2[i];
return arr;
}


static bool FluxGreater(const SourceData& x, const SourceData& y) { return x.flux>y.flux; }

void SourceDataArray::SortByFlux()
{
sort(m_dataArr.begin(), m_dataArr.end(), FluxGreater);
}


static SourceData errData;

int SourceDataArray::IndexOf(string label) const
{
int count = Count();
for (int i=0; i<count; ++i)
	if (m_dataArr[i].label==label)
		return i;
return -1;
}

const SourceData& SourceDataArray::operator[](string label) const
{
int index = IndexOf(label);
return index<0 ? errData : m_dataArr[index];
}

SourceData& SourceDataArray::operator[](string label)
{
int index = IndexOf(label);
return index<0 ? errData : m_dataArr[index];
}

void SourceDataArray::Print(ostream& oFile) const
{
int count = Count();
for (int i=0; i<count; ++i)
	m_dataArr[i].Print(oFile);
}


// Append the default value (zero) if the string is not all spaces or empty
// Return false otherwise
static bool AppendDefault(char* buff, stringstream& str)
{
	while (*buff && *buff<=32) {
		++buff;
	}
	if (!*buff)
		return false;
	strcat(buff, " 0.0 0 0 0");
	str.str(buff);
	return true;
}

extern SourceDataArray ReadSourceFile(const char* srcfilename)
{
SourceDataArray srcArr;
if (!strcmp(srcfilename, c_missingFileName))
	return srcArr;
ifstream infile(srcfilename);
if (!infile.is_open())
	return srcArr;
char buffer[1024];
while (!infile.eof()) {
	infile.getline(buffer, 1024);
	if (buffer[0]!='!') {
		stringstream str(ios_base::in);
		cout << "## " << buffer << endl;
		if (AppendDefault(buffer, str)) {
			SourceData srcData;
			str >> srcData.flux >> srcData.srcL >> srcData.srcB >> srcData.index
			>> srcData.fixflag >> srcData.minTS >> srcData.label >> srcData.loclimit >> srcData.typefun >> srcData.par2 >> srcData.par3;
			srcData.minTS = srcData.minTS * srcData.minTS;
			srcArr.Append(srcData);
			}
		}
	}
infile.close();
return srcArr;
}



/// ////////////////////////////////////////////
///
/// MapList and MapData classes
///
///


// Remove leading and trailing spaces
static string Trim(string s)
{
int len = s.length();
if (!len)
	return s;
int a = 0;
while (a<len && s.at(a)<=32)
	++a;
if (a==len)
	return s="";
int z = len-1;
while (z>0 && s.at(z)<=32)
	--z;
return s.substr(a, z-a+1);
}




bool MapList::Read(const char* listFileName)
{
m_count = 0;
m_ctslist.clear();
m_explist.clear();
m_gaslist.clear();
m_thetaArr.clear();
if (!strcmp(listFileName, c_missingFileName))
	return true;

string line;
string ctsName;
string expName;
string gasName;
double theta;
double galCoeff;
double isoCoeff;
ifstream ifile(listFileName);
if (ifile.is_open()) {
	while (!ifile.eof()) {
		getline(ifile, line);
		line = Trim(line);
		if (line.length() && line.at(0)!='!') {
			++m_count;
			stringstream str(ios_base::in);
			str.str(line.c_str());
			str >> ctsName >> expName >> gasName >> theta >> galCoeff >> isoCoeff;
			m_ctslist.push_back(ctsName);
			m_explist.push_back(expName);
			m_gaslist.push_back(gasName);
			m_thetaArr.push_back(theta);
			m_galCoeffArr.push_back(galCoeff);
			m_isoCoeffArr.push_back(isoCoeff);
			}
		}
	ifile.close();
	return true;
	}
cerr << "Error opening " << listFileName << endl;
return false;
}



/// //////////////////////////////////////////////////////////////////////////////////////




void MapCoeff::Clear()
{
m_length = 0;
m_thetaArr = 0;
m_galCoeffs = 0;
m_isoCoeffs = 0;
}

void MapCoeff::Destroy()
{
delete[] m_thetaArr;
delete[] m_galCoeffs;
delete[] m_isoCoeffs;
}


void MapCoeff::Copy(const MapCoeff& another)
{
Clear(); /// Clear the arrays in case m_length==0
m_length = another.m_length;
if (m_length) {
	m_thetaArr = new double[m_length];
	m_galCoeffs = new double[m_length];
	m_isoCoeffs = new double[m_length];
	for (int i=0; i<m_length; ++i) {
		m_thetaArr[i] = another.m_thetaArr[i];
		m_galCoeffs[i] = another.m_galCoeffs[i];
		m_isoCoeffs[i] = another.m_isoCoeffs[i];
		}
	}
}


bool MapCoeff::Load(const MapList& maplist)
{
Destroy();
Clear();
m_length = maplist.Count();
if (!m_length)
	return false;

m_thetaArr = new double[m_length];
m_galCoeffs = new double[m_length];
m_isoCoeffs = new double[m_length];
for (int i=0; i<m_length; ++i) {
	m_thetaArr[i] = maplist.Theta(i);
	m_galCoeffs[i] = maplist.GalCoeff(i);
	m_isoCoeffs[i] = maplist.IsoCoeff(i);
	}
return true;
}


MapCoeff::MapCoeff(double theta, double galCoeff, double isoCoeff)
{
m_length = 1;
m_thetaArr = new double[1];
m_galCoeffs = new double[1];
m_isoCoeffs = new double[1];
m_thetaArr[0] = theta;
m_galCoeffs[0] = galCoeff;
m_isoCoeffs[0] = isoCoeff;
}



/// //////////////////////////////////////////////////////////////////////////////////////



MapMaps::MapMaps(const AgileMap& ctsMap, const AgileMap& expMap, const AgileMap& gasMap)
{
m_mapCount = 1;
m_ctsMaps = new AgileMap[m_mapCount];
m_expMaps = new AgileMap[m_mapCount];
m_gasMaps = new AgileMap[m_mapCount];
m_ctsMaps[0] = ctsMap;
m_expMaps[0] = expMap;
m_gasMaps[0] = gasMap;
m_alignedMaps = true;
m_loadFlags = new int[m_mapCount];
m_countsArr = new int[m_mapCount];
m_energyInf = expMap.GetEmin();
m_energySup = expMap.GetEmax();
double counts = 0;
for (int r=0; r<ctsMap.Rows(); ++r)
	for (int c=0; c<ctsMap.Cols(); ++c)
		counts += ctsMap(r, c);
m_countsArr[0] = int(counts);
m_loadFlags[0] = ctsMap.SameParamsAs(expMap) ? 0 : 16;
}


void MapMaps::Clear()
{
m_mapCount = 0;
m_ctsMaps = 0;
m_expMaps = 0;
m_gasMaps = 0;
m_loadFlags = 0;
m_countsArr = 0;
m_alignedMaps = true;
m_energyInf = 0;
m_energySup = 0;
}

void MapMaps::Destroy()
{
delete[] m_ctsMaps;
delete[] m_expMaps;
delete[] m_gasMaps;
delete[] m_countsArr;
delete[] m_loadFlags;
}


void MapMaps::Copy(const MapMaps& another)
{
Clear(); /// Clear the arrays in case m_mapCount==0
m_mapCount = another.m_mapCount;
m_alignedMaps = another.m_alignedMaps;
m_energyInf = another.m_energyInf;
m_energySup = another.m_energySup;
if (m_mapCount) {
	m_ctsMaps = new AgileMap[m_mapCount];
	m_expMaps = new AgileMap[m_mapCount];
	m_gasMaps = new AgileMap[m_mapCount];
	m_countsArr = new int[m_mapCount];
	m_loadFlags = new int[m_mapCount];
	for (int i=0; i<m_mapCount; ++i) {
		m_ctsMaps[i] = another.m_ctsMaps[i];
		m_expMaps[i] = another.m_expMaps[i];
		m_gasMaps[i] = another.m_gasMaps[i];
		m_countsArr[i] = another.m_countsArr[i];
		m_loadFlags[i] = another.m_loadFlags[i];
		}
	}
}





bool MapMaps::SetCtsMap(const AgileMap& ctsMap, int i)
{
if (i<0 || i>=m_mapCount) {
	cerr << "MapMaps::ERROR: i=" << i << endl;
	return false;
	}
if (!ctsMap.SameParamsAs(m_expMaps[i])) {
	cerr << "MapMaps::ERROR: !SameParamsAs(" << i << ")" << endl;
	m_loadFlags[i] = m_loadFlags[i] | 16;
	return false;
	}
else
	m_loadFlags[i] = m_loadFlags[i] & ~16;
m_loadFlags[i] = m_loadFlags[i] | 1;
if (m_alignedMaps && m_mapCount>1)
	m_alignedMaps = ctsMap.AlignedTo(m_ctsMaps[(i+1)%m_mapCount]);

double counts = 0;
for (int r=0; r<ctsMap.Rows(); ++r)
	for (int c=0; c<ctsMap.Cols(); ++c)
		counts += ctsMap(r, c);
m_countsArr[i] = int(counts);
m_ctsMaps[i] = ctsMap;
return true;
}




bool MapMaps::ReplaceCtsMaps(AgileMap* ctsMapArr)
{
for (int m=0; m<m_mapCount; ++m) {
	m_loadFlags[m] = m_loadFlags[m] | 1;
	const AgileMap& ctsMap = ctsMapArr[m];
	if (!ctsMap.SameParamsAs(m_expMaps[m])) {
		cerr << "MapMaps::ERROR: !SameParamsAs(" << m << ")" << endl;
		m_loadFlags[m] = m_loadFlags[m] | 16;
		return false;
		}
	else
		m_loadFlags[m] = m_loadFlags[m] & ~16;

	if (m_alignedMaps && m)
		m_alignedMaps = ctsMap.AlignedTo(m_ctsMaps[m-1]);
	double counts = 0;
	for (int r=0; r<ctsMap.Rows(); ++r)
		for (int c=0; c<ctsMap.Cols(); ++c)
			counts += ctsMap(r, c);
	m_countsArr[m] = int(counts);
	}
//delete[] m_ctsMaps; //AB 2015-11-01
*m_ctsMaps = *ctsMapArr;
return true;
}




bool MapMaps::Load(const MapList& maplist, bool skipCts)
{
Destroy();
Clear();
m_mapCount = maplist.Count();
if (!m_mapCount)
	return false;

m_ctsMaps = new AgileMap[m_mapCount];
m_expMaps = new AgileMap[m_mapCount];
m_gasMaps = new AgileMap[m_mapCount];
m_countsArr = new int[m_mapCount];
m_loadFlags = new int[m_mapCount];

m_energyInf = 0;
m_energySup = 0;
m_alignedMaps = true;

bool mapError = false;
for (int i=0; i<m_mapCount && !mapError; ++i) {
	AgileMap& ctsMap = m_ctsMaps[i];
	AgileMap& expMap = m_expMaps[i];
	AgileMap& gasMap = m_gasMaps[i];
	const char* ctsName = maplist.CtsName(i);
	const char* expName = maplist.ExpName(i);
	const char* gasName = maplist.GasName(i);
	bool ctsErr = false;
	bool expErr = false;
	bool gasErr = false;
	m_loadFlags[i] = 0;

	expErr = expMap.Read(expName) || !strcmp(gasName, "None");
	if (expErr)
		cerr << "Errors reading the exposure map" << endl;

	m_countsArr[i] = 0;
	if (skipCts) {
		ctsMap = expMap;
		ctsMap.Zero();
		m_loadFlags[i] = m_loadFlags[i] | 1;
		}
	else {
		ctsErr = ctsMap.Read(ctsName) || !strcmp(ctsName, "None");
		if (ctsErr)
			cerr << "Errors reading the counts map" << endl;
		else {
			double counts = 0;
			for (int r=0; r<ctsMap.Rows(); ++r)
				for (int c=0; c<ctsMap.Cols(); ++c)
					counts += ctsMap(r, c);
			m_countsArr[i] = int(counts);
			}
		}

	gasErr = gasMap.Read(gasName) || !strcmp(gasName, "None");
	if (gasErr)
		cerr << "Errors reading the gas map" << endl;

	m_loadFlags[i] = skipCts | (ctsErr<<1) | (expErr<<2) | (gasErr<<3);

	mapError = mapError || ctsErr || expErr || gasErr;
	if (!mapError) {
		if (i==0 || m_energyInf>expMap.GetEmin())
			m_energyInf = expMap.GetEmin();
		if (i==0 || m_energySup<expMap.GetEmax())
			m_energySup = expMap.GetEmax();
		if (!skipCts && !ctsMap.SameParamsAs(expMap)) {
			m_loadFlags[i] = m_loadFlags[i] | 16;
			mapError = true;
			cerr << "Params mismatch between " << ctsName << " and " << expName << endl;
			}
		if (!expMap.AlignedTo(gasMap)) {
			m_loadFlags[i] = m_loadFlags[i] | 32;
			mapError = true;
			cerr << "Alignment mismatch between " << expName << " and " << gasName << endl;
			}
		if (m_alignedMaps && i>0)
			m_alignedMaps = expMap.AlignedTo(m_expMaps[0]);
		}
	}
cout << endl;
if (mapError) {
	cerr << "Errors reading the map list" << endl;
	Destroy();
	Clear();
	}
else if (m_mapCount==1)
	cout << "Using a single map set" << endl;
else if (m_alignedMaps)
	cout << "Using " << m_mapCount << " aligned map sets" << endl;
else
	cout << "Using " << m_mapCount << " NON aligned map sets" << endl;
return !mapError;
}



bool MapData::Load(const MapList& maplist, bool skipCts)
{
bool ok = MapMaps::Load(maplist, skipCts);
ok = MapCoeff::Load(maplist) && ok;
return ok && m_length==m_mapCount;
}



void MapData::Print() const
{
int count = m_length;
if (count>m_mapCount)
	count = m_mapCount;
cout << "Maps cts exp gas theta eMin eMax galCoeff isoCoeff" << endl;
for (int i=0; i<count; ++i) {
	const char* ctsName;
	const char* expName;
	const char* gasName;
	const char* ctsXexp = " ";
	const char* expXgas = " ";
	int flags = m_loadFlags[i];
	if (flags&2)
		ctsName = "Error";
	else if (flags&1)
		ctsName = "None";
	else
		ctsName = m_ctsMaps[i].GetFileName();
	if (flags&4)
		expName = "Error";
	else
		expName = m_expMaps[i].GetFileName();
	if (flags&8)
		gasName = "Error";
	else
		gasName = m_gasMaps[i].GetFileName();
	if (flags&16)
		ctsXexp = "<X>";
	if (flags&32)
		expXgas = "<X>";
	cout << i+1 << " " << ctsName << ctsXexp << expName << expXgas << gasName << " " << Theta(i) << " " << GalCoeff(i) << " " << IsoCoeff(i) << endl;
	}
if (m_mapCount>1) {
	if (m_alignedMaps)
		cout << "All maps aligned" << endl;
	else
		cout << "Maps not aligned" << endl;
	}
cout << "Energy from " << m_energyInf << " to " << m_energySup  << endl;
}




/// ////////////////////////////////////////////
///
/// ExtList and ExtData classes
///
///


bool ExtList::Read(const char* listFileName)
{
Clear();
if (!strcmp(listFileName, c_missingFileName))
	return true;

bool syntaxOK = true;
int lineNum = 0;
string line;
ifstream ifile(listFileName);
if (ifile.is_open()) do {
	getline(ifile, line);
	line = Trim(line);
	++lineNum;
	if (line.length()!=0 && line.at(0)!='!') {
		stringstream str(ios_base::in);
		str.str(line.c_str());
		string name;
		double flux;
		double index;
		str >> name >> flux >> index;
		m_extNames.push_back(name);
		m_fluxArr.push_back(flux);
		m_indexArr.push_back(index);
		++m_extCount;
		int mapCount = 0;
		while (!str.eof()) {
			string map;
			str >> map;
			++mapCount;
			m_extMaps.push_back(map);
			}
		if (m_extCount==1)
			m_mapCount = mapCount;
		else if (m_mapCount!=mapCount) {
			syntaxOK = false;
			cerr << "Error in " << listFileName << " line " << lineNum << ": " << m_mapCount << " expected, " << mapCount << " found" << endl;
			}
		}
	if (!syntaxOK) {
		Clear();
		ifile.close();
		break;
		}
	} while (!ifile.eof());
else {
	syntaxOK = false;
	cerr << "Error opening file " << listFileName << endl;
	}
return syntaxOK;
}





bool ExtData::Load(const ExtList& extlist)
{
Destroy();
Clear();
m_mapCount = extlist.MapCount();
m_extCount = extlist.ExtCount();
cout << "Reading " << m_mapCount*m_extCount << " map for " << m_extCount << " extended source";
if (m_mapCount>1)
	cout << "s";
cout << endl;
if (!m_mapCount || !m_extCount)
	return true;

m_nameOfs = new int[m_extCount];
m_nameBuffSize = 0;
for (int i=0; i<m_extCount; ++i) {
	m_nameOfs[i] = m_nameBuffSize;
	m_nameBuffSize += strlen(extlist.ExtName(i))+1;
	}
m_nameBuff = new char[m_nameBuffSize];
for (int i=0; i<m_extCount; ++i)
	strcpy(&m_nameBuff[m_nameOfs[i]], extlist.ExtName(i));

m_extArr = new AgileMap[m_mapCount*m_extCount];
m_fluxArr = new double[m_extCount];
m_indexArr = new double[m_extCount];
bool errors = false;
for (int m=0; m<m_mapCount && !errors; ++m)
	for (int i=0; i<m_extCount && !errors; ++i) {
		/// cerr << "m=" << m << ", i=" << i << ", Idx(m, i)=" << Idx(m, i) << endl;
		errors = m_extArr[Idx(m, i)].Read(extlist.ExtFileName(m, i));
		}
for (int i=0; i<m_extCount; ++i) {
	m_fluxArr[i] = extlist.ExtFlux(i);
	m_indexArr[i] = extlist.ExtIndex(i);
	}
return !errors;
}

void ExtData::Copy(const ExtData& another)
{
m_mapCount = another.m_mapCount;
m_extCount = another.m_extCount;
if (!m_mapCount || !m_extCount)
	return;
m_extArr = new AgileMap[m_mapCount*m_extCount];
m_fluxArr = new double[m_extCount];
m_indexArr = new double[m_extCount];

m_nameBuffSize = another.m_nameBuffSize;
m_nameBuff = new char[m_nameBuffSize];
m_nameOfs = new int[m_extCount];
memcpy(m_nameBuff, another.m_nameBuff, m_nameBuffSize);
memcpy(m_nameOfs, another.m_nameOfs, m_extCount*sizeof(int));

for (int i=0; i<m_mapCount*m_extCount; ++i)
	m_extArr[i] = another.m_extArr[i];
for (int i=0; i<m_extCount; ++i) {
	m_fluxArr[i] = another.m_fluxArr[i];
	m_indexArr[i] = another.m_indexArr[i];
	}
}
