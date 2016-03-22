//
// C++ Implementation: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <TGraph.h>
#include <TMatrixDEigen.h>
#include <TRandom3.h>
#ifndef ROOT_TROOT
#include "TROOT.h"
#endif

#include <iomanip>
#include <iostream>
#include <fstream>
#include <ctime>

#include "wcstrig.h"

#include "RoiMulti5.h"

using namespace std;


const int FILENAME_SIZE = 1025;

Double_t EvalExposure(Double_t srcl, Double_t srcb, const AgileMap& expmap)
{
int srcR, srcC;
bool inside = expmap.GetRowCol(srcl, srcb, &srcR, &srcC);
if (inside)
	return expmap(srcR, srcC) / expmap.Area(srcR, srcC);
return 0;
}



void AlikeDiffMap::Copy(const AlikeDiffMap &another, bool all)
{
if (all) {
	AgileMap::operator=(another);
	AlikeNorm::operator=(another);
	}
m_index = another.m_index;
m_coeff = another.m_coeff;
m_coefferr = another.m_coefferr;
m_coefflo = another.m_coefflo;
m_coeffhi = another.m_coeffhi;
m_fixflag = another.m_fixflag;
}



void AlikeDiffMap::BuildFrom(const AgileMap& map, double eInf, double eSup, double index)
{
AgileMap::operator=(map);
m_index = index;
m_coeff = 1;
m_coefferr = 0;
m_coefflo = 0;
m_coeffhi = 0;
m_fixflag = 1;
SetEnergyRange(eInf, eSup);
UpdateNorm();
}


void AlikeDiffMap::Expose(const AgileMap& expMap)
{
MatD& gas = *this;
int size = expMap.Size();
for (int i=0; i<size; ++i)
	gas[i] *= expMap[i];
}



/**
void AlikeDiffMap::UpdNorm()
{
m_normFactor = (pow(GetEmin(), 1.0-m_index)-pow(GetEmax(), 1.0-m_index)) / pow(100.0, 1.0-m_index);
}
*/


void AlikeExtMap::Copy(const AlikeExtMap &another, bool all)
{
if (all)
	AlikeDiffMap::Copy(another, all);
m_TS=another.m_TS;
m_baryl=another.m_baryl;
m_baryb=another.m_baryb;
}

void AlikeExtMap::EvalBarycenter()
{
int rows = Rows();
int cols = Cols();
const AlikeExtMap& map = *this;
double sum = 0;
double isum = 0;
double jsum = 0;
for (int i=0; i<rows; ++i)
	for (int j=0; j<cols; ++j) {
		double val = map(i, j);
		isum += val*i;
		jsum += val*j;
		sum += val;
		}
int baryi = int(isum/sum+0.5);
int baryj = int(jsum/sum+0.5);
m_baryl = l(baryi, baryj);
m_baryb = b(baryi, baryj);
cerr << "Barycenter: " << m_baryl << ", " << m_baryb << endl;
}

void AlikeExtMap::SetTS(double ts)
{
	m_TS = ts<0 ? 0 : ( gROOT->GetVersionInt()<53000 ? ts : 2.0*ts);
}

void AlikeExtMap::PrintSqrtTS() const
{
cout << "sqrt(TS) = " << sqrt(GetTS()) << endl;
}

/// /////////////////////////////////////////////////////////
///
/// Class AlikeSourceMap
///



void AlikeSourceMap::WritePolygon(string fileprefix) const
{
if (m_polygon.Sides()) {
	string fileName(fileprefix);
	fileName += ".con";
	ofstream conout(fileName.c_str());
	/// Printing the polygon vertices (the first and the last being the same)
	for (int i=0; i<=m_polygon.Sides(); ++i)
		conout << m_polygon[i].x << " " << m_polygon[i].y << endl;
	conout.close();
	}
}

void AlikeSourceMap::WriteEllipse(string fileprefix) const
{
const double numPoints = 40;

if (m_ellipse.horAxis>0 && m_ellipse.verAxis>0) {
	string fileName(fileprefix);
	fileName += ".ellipse.con";
	ofstream conout(fileName.c_str());
	double aQ = m_ellipse.horAxis*m_ellipse.horAxis;
	double bQ = m_ellipse.verAxis*m_ellipse.verAxis;
	double cosb = cos(m_ellipse.center.y*DEG2RAD);
	double stepAngle = 2*M_PI/numPoints;
	for (int i=0; i<=numPoints; ++i) {
		double angle = i*stepAngle-m_ellipse.attitude;
		double sint = sin(angle);
		double sinQt = sint*sint;
		double cost = cos(angle);
		double cosQt = cost*cost;
		double r = 1.0/sqrt(cosQt/aQ+sinQt/bQ);
		double x = r*cos(i*stepAngle)/cosb;
		double y = r*sin(i*stepAngle);
		conout << m_ellipse.center.x+x << " " << m_ellipse.center.y+y << endl;
		}
	fileName = fileprefix + ".reg";
	ofstream conout2(fileName.c_str());
	conout2 << "galactic" << endl;
	conout2 << "ellipse(" << m_ellipse.center.x << "," << m_ellipse.center.y << ","
			<< m_ellipse.horAxis << "," << m_ellipse.verAxis << "," << -m_ellipse.attitude*RAD2DEG << ")" << endl;
	}
}

void AlikeSourceMap::SetTS(Double_t inTS)
{
	m_ts = inTS<0 ? 0 : ( gROOT->GetVersionInt()<53000 ? inTS : 2.0 * inTS);
}

void AlikeSourceMap::PrintSqrtTS() const
{
cout << "sqrt(TS) = " << sqrt(GetTS()) << endl;
}

void AlikeSourceMap::SetULCL(Double_t ulcl)
{
	m_ulcl = ulcl;
}

void AlikeSourceMap::SetLocCL(Double_t loccl)
{
	m_loccl = loccl;
}

/// ////////////////////////////////////////////////////
///
/// class FitParameters
///
/// Class for remapping the model parameters indices
///



void FitParameters::Init(int count, bool allVariables)
{
delete[] m_params;
delete[] m_vars;
delete[] m_map;
m_count = count;
m_params = new double[m_count];
m_map = new int[m_count];
m_vars = new bool[m_count];
m_mapCount = allVariables ? m_count : 0;
for (int i=0; i<m_count; ++i) {
	m_params[i] = 0;
	m_map[i] = 0;
	m_vars[i] = allVariables;
	}
}


void FitParameters::SetValue(int i, double value, bool makeVariable)
{
m_params[i] = value;
if (m_vars[i]!=makeVariable) {
	m_vars[i] = makeVariable;
	UpdateMap();
	}
}


void FitParameters::UpdateValues(const double* iPars)
{
for (int i=0; i<m_mapCount; ++i)
	m_params[m_map[i]] = iPars[i];
}


void FitParameters::FillParArray(const double* varVals, double* allVals) const
{
for (int i=0; i<m_count; ++i)
	allVals[i] = m_params[i];
for (int i=0; i<m_mapCount; ++i)
	allVals[m_map[i]] = varVals[i];
}


void FitParameters::SetModelValues(TF1& model) const
{
int varIndex = 0;
for (int i=0; i<m_count; ++i)
	if (m_vars[i])
		model.SetParameter(varIndex++, m_params[i]);
}

void FitParameters::GetModelValues(const TF1& model)
{
for (int i=0; i<m_mapCount; ++i)
	m_params[m_map[i]] = model.GetParameter(i);
}


void FitParameters::UpdateMap()
{
int varIndex = 0;
for (int i=0; i<m_count; ++i)
	if (m_vars[i])
		m_map[varIndex++] = i;
}




RoiMulti* RoiMulti::sm_singleton = 0;

int RoiMulti::s_fitterIterations = 0;

unsigned* RoiMulti::s_fitterChanges = 0;


RoiMulti::RoiMulti():

	m_countsHist("Counts","Counts",0,0,0), m_model(),
	m_ranal(0.0), m_ulcl(0.0), m_loccl(0.0), m_lCenter(0.0), m_bCenter(0.0),


	m_psfTab(),
	m_mapCount(0),
	m_ctsMaps(0),
	m_expMaps(0),
	m_gasMaps(0),
	m_alignedMaps(false),

	m_counts(0),
	m_thetaArr(0),
	m_galCoeffs(0),
	m_isoCoeffs(0),

	m_energyInf(0.0),
	m_energySup(0.0),

	m_inSrcDataArr(0),
	m_outSrcDataArr(0),
	
	m_galMode(DiffFixed),
	m_isoMode(DiffFixed),
	m_galSrc(0),
	m_isoSrc(0),

	m_galParCount(0),
	m_isoParCount(0),
	m_diffParCount(0),

	m_extData(),
	m_extCount(0),
	m_extSrcArr(0),
	m_extFlux(0),
	m_extIndex(0),

	m_srcCount(0),
	m_sources(0),
	m_srcLimits(0),
	m_fitInfo(0),

	m_sourceParOffset(0),
	m_logFile(0),

	m_binList(0),
	m_binAddr(0),

	m_binCount(0),
	m_fluxScaleFactor(0.0),
	m_fluxScaleFactorMulInv(0.0),
	m_fluxUpperBound(0.0)
{
	m_countsHist.Sumw2();
if (sm_singleton)
	throw "RoiMulti allocation Error";
else
	sm_singleton = this;
}


RoiMulti::~RoiMulti()
{
CleanMaps();
CleanSources();
delete[] m_binList;
delete[] m_binAddr;
CloseLogFile();
sm_singleton = 0;
}



void RoiMulti::CleanMaps()
{
delete[] m_ctsMaps;
delete[] m_expMaps;
delete[] m_gasMaps;
delete[] m_counts;
delete[] m_thetaArr;
delete[] m_galCoeffs;
delete[] m_isoCoeffs;

delete[] m_extSrcArr;
delete[] m_extFlux;
delete[] m_extIndex;


m_ctsMaps = 0;
m_expMaps = 0;
m_gasMaps = 0;
m_extSrcArr = 0;
m_extFlux = m_extIndex = 0;

m_extCount = 0;
m_counts = 0;
m_galCoeffs = m_isoCoeffs = m_thetaArr = 0;
m_mapCount = 0;
m_diffParCount = 0;
m_sourceParOffset = 0;
}

void RoiMulti::CleanSources()
{
delete[] m_galSrc;
delete[] m_isoSrc;
m_galSrc = m_isoSrc = 0;

delete[] m_sources;
delete[] m_srcLimits;
delete[] m_fitInfo;
m_sources = 0;
m_srcLimits = 0;
m_fitInfo = 0;
m_srcCount = 0;
}




bool RoiMulti::SetPsf(const char* psfFileName, const char* raeffFileName, const char* edpFileName)
{
CleanSources();
CleanMaps();
return m_psfTab.Read(psfFileName, raeffFileName, edpFileName);
}




/// /////////////////////////////////////////////
///
/// Simulation
///
///


class PoissonStat
{
public:
	PoissonStat(); ///  { m_random.SetSeed(0); }
	~PoissonStat() {}
	int Poisson(double mean) { return m_random.Poisson(mean); }
	unsigned GetSeed() { return m_random.GetSeed(); }
	void SetSeed(unsigned seed) { m_random.SetSeed(seed); }
private:
	TRandom3 m_random;
};

PoissonStat::PoissonStat()
{
time_t seed = time(0);
m_random.SetSeed(seed);
}

static PoissonStat s_poissonStat;

extern unsigned GetSeed()
{
return s_poissonStat.GetSeed();
}

extern unsigned SetSeed(unsigned seed)
{
if (!seed)
	seed = time(0);
unsigned oldSeed = s_poissonStat.GetSeed();
s_poissonStat.SetSeed(seed);
return oldSeed;
}



static void AddMap(AgileMap& target, double factor, const AgileMap& add)
{
int rows = target.Rows();
int cols = target.Cols();
if (factor==0 || rows!=add.Rows() || cols!=add.Cols())
	return;
for (int i=0; i<rows; ++i)
	for (int j=0; j<cols; ++j)
		target(i,j) += factor*add(i,j);
}




AgileMap* RoiMulti::NewModelArray(const SourceDataArray& srcDataArr)
{
if (!m_mapCount) {
	cerr << "No maps loaded: Simulation not possible" << endl;
	return 0;
	}
SetSources(srcDataArr, 10, 2, 5.99);

cout << "Generating the model..." << flush;
AgileMap* modelArr = new AgileMap[m_mapCount];
for (int map=0; map<m_mapCount; ++map) {
	AgileMap model(m_isoSrc[map]);
	model *= m_isoCoeffs[map];
	if (m_galCoeffs[map] != 0)
		model += m_galSrc[map]*m_galCoeffs[map];
	for (int s=0; s<m_srcCount; ++s) {
		AlikeSourceMap& source = m_sources[map*m_srcCount+s];
		/// Add(model, source.GetFlux()*source.GetExp()*source.GetNormFactor(), source);
		AddMap(model, source.GetFlux()*source.GetExp()*source.GetNormFactor(), source);
		}
	modelArr[map] = model;
	}
cout << "Done" << endl;
return modelArr;
}



AgileMap* RoiMulti::NewSimulationArray(const SourceDataArray& srcDataArr)
{
AgileMap* modelArr = NewModelArray(srcDataArr);
if (!modelArr)
	return 0;

/// gRandom->SetSeed(0);

AgileMap* simArr = new AgileMap[m_mapCount];
for (int map=0; map<m_mapCount; ++map) {
	const AgileMap& model = modelArr[map];
	AgileMap& sim = simArr[map];
	sim = model;
	int rows = sim.Rows();
	int cols = sim.Cols();
	for (int i=0; i<rows; ++i)
		for (int j=0; j<cols; ++j)
			sim(i,j) = s_poissonStat.Poisson(model(i,j));
	}
delete[] modelArr;
return simArr;
}




bool RoiMulti::WriteSimulatedCounts(const SourceDataArray& srcDataArr, const char* outfilename)
{
AgileMap* simArr = NewSimulationArray(srcDataArr);
if (!simArr)
	return false;
char name[FILENAME_SIZE];
bool errors = false;
for (int map=0; map<m_mapCount; ++map) {
	const char* format = "%s_%d.sim.cts.gz";
	if (m_mapCount==1)
		format = "%s.sim.cts.gz";
	sprintf(name, format, outfilename, map+1);
	if (simArr[map].Write(name)) {
		cerr << "Error writing simulated counts map " << name << endl;
		errors = true;
		}
	}
delete[] simArr;
return !errors;
}


bool RoiMulti::WriteSimulatedCounts(const SourceDataArray& srcDataArr, const MapList& fileNames)
{
AgileMap* simArr = NewSimulationArray(srcDataArr);
if (!simArr)
	return false;
bool errors = false;
for (int map=0; map<m_mapCount; ++map) {
	if (simArr[map].Write(fileNames.CtsName(map))) {
		cerr << "Error writing simulated counts map " << fileNames.CtsName(map) << endl;
		errors = true;
		}
	}
delete[] simArr;
return !errors;
}


bool RoiMulti::WriteModelCounts(const SourceDataArray& srcDataArr, const char* outfilename)
{
AgileMap* modelArr = NewModelArray(srcDataArr);
if (!modelArr)
	return false;
char name[FILENAME_SIZE];
bool errors = false;
for (int map=0; map<m_mapCount; ++map) {
	const char* format = "%s_%d.model.cts.gz";
	if (m_mapCount==1)
		format = "%s.model.cts.gz";
	sprintf(name, format, outfilename, map+1);
	if (modelArr[map].Write(name)) {
		cerr << "Error writing model counts map " << name << endl;
		errors = true;
		}
	}
delete[] modelArr;
return !errors;
}


bool RoiMulti::WriteModelCounts(const SourceDataArray& srcDataArr, const MapList& fileNames)
{
AgileMap* modelArr = NewModelArray(srcDataArr);
if (!modelArr)
	return false;
bool errors = false;
for (int map=0; map<m_mapCount; ++map) {
	char modelname[FILENAME_SIZE] = "model_";
	strcat(modelname,fileNames.CtsName(map));
	if (modelArr[map].Write(modelname)) {
		cerr << "Error writing model counts map " << modelname << endl;
		errors = true;
		}
	}
delete[] modelArr;
return !errors;
}




static bool CheckDiffMode(const char* paramName, DiffMode paramValue)
{
if (paramValue>=DiffFixed && paramValue<=DiffTogether)
	return true;
cerr << paramName << " must range between " << DiffFixed << " and " << DiffTogether << endl;
return false;
}


bool RoiMulti::SetMaps(const MapData& mapData, DiffMode galMode, DiffMode isoMode)
{
CleanSources();
CleanMaps();

m_mapCount = mapData.Count();
if (!m_mapCount) {
	cerr << "No Maps found" << endl;
	return false;
	}
if (!(CheckDiffMode("galmode", galMode) && CheckDiffMode("isomode", isoMode)))
	return false;

m_galMode = galMode;
m_isoMode = isoMode;

m_alignedMaps = mapData.Aligned();
m_energyInf = mapData.EnergyInf();
m_energySup = mapData.EnergySup();

m_ctsMaps = new AgileMap[m_mapCount];
m_expMaps = new AgileMap[m_mapCount];
m_gasMaps = new AgileMap[m_mapCount];
m_thetaArr = new double[m_mapCount];
m_galCoeffs = new double[m_mapCount];
m_isoCoeffs = new double[m_mapCount];
m_counts = new int[m_mapCount];

for (int map=0; map<m_mapCount; ++map) {
	m_thetaArr[map] = mapData.Theta(map);
	m_galCoeffs[map] = mapData.GalCoeff(map);
	m_isoCoeffs[map] = mapData.IsoCoeff(map);
	m_counts[map] = mapData.MapCounts(map);

	m_ctsMaps[map] = mapData.CtsMap(map);

	/// zzz temporary
	double totCts = 0;
	for (int r=0; r<m_ctsMaps[map].Rows(); ++r)
		for (int c=0; c<m_ctsMaps[map].Cols(); ++c)
			totCts += m_ctsMaps[map](r,c);
	cout << "Total CTS = " << totCts << endl;


	m_expMaps[map] = mapData.ExpMap(map);
	m_gasMaps[map] = mapData.GasMap(map);
	}

m_galSrc = new AlikeDiffMap[m_mapCount];
m_isoSrc = new AlikeDiffMap[m_mapCount];
for (int map=0; map<m_mapCount; ++map) {
	AlikeDiffMap& galMap = m_galSrc[map];
	galMap.BuildFrom(m_gasMaps[map], m_energyInf, m_energySup);	/// Default spectral index

	/// galMap.TMatrixD::operator=(ElementMult(galMap, m_expMaps[map]));
	galMap.Expose(m_expMaps[map]);

	/// zzz temporary
	double totExp = 0;
	for (int r=0; r<galMap.Rows(); ++r)
		for (int c=0; c<galMap.Cols(); ++c)
			totExp += galMap(r,c);
	cout << "Total exposure = " << totExp << endl;

	FixFlag fixFlag = FluxFree;
	if (m_galMode==DiffFixed || (m_galMode==DiffDefault && m_galCoeffs[map]>0))
		fixFlag = AllFixed;
	galMap.SetFixflag(fixFlag);
	m_galCoeffs[map] = fabs(m_galCoeffs[map]);
	galMap.SetCoeff(m_galCoeffs[map]);
	cout << "Gasmap " << map+1 << " " << galMap.GetFileName() << " index = " << galMap.GetIndex() << endl;

	AlikeDiffMap& isoMap = m_isoSrc[map];
	isoMap.BuildFrom(m_expMaps[map], m_energyInf, m_energySup);	/// Default spectral index

	isoMap *= 1.0e-5;

	fixFlag = FluxFree;
	if (m_isoMode==DiffFixed || (m_isoMode==DiffDefault && m_isoCoeffs[map]>0))
		fixFlag = AllFixed;
	isoMap.SetFixflag(fixFlag);
	m_isoCoeffs[map] = fabs(m_isoCoeffs[map]);
	isoMap.SetCoeff(m_isoCoeffs[map]);
	cout << "Isomap " << map+1 << " " << isoMap.GetFileName() << " index = " << isoMap.GetIndex() << endl;

	/// zzz temporary
	double totIso = 0;
	for (int r=0; r<isoMap.Rows(); ++r)
		for (int c=0; c<isoMap.Cols(); ++c)
			totIso += isoMap(r,c);
	cout << "Total GAL = " << totIso << endl;

	}
return true;
}





bool RoiMulti::SetExtendedSources(const ExtData& extData)
{
m_extData = extData;
if (!m_extData.ExtCount()) {
	cout << "No extended sources considered" << endl;
	return true;
	}
if (m_extData.MapCount()!=m_mapCount) {
	cerr << "ERROR: The number of counts maps and extended sources maps differ" << endl;
	return false;
	}
delete[] m_extSrcArr;
delete[] m_extFlux;
delete[] m_extIndex;
m_extCount = m_extData.ExtCount();
if (m_extCount) {
	m_extSrcArr = new AlikeExtMap[m_mapCount*m_extCount];
	m_extFlux = new double[m_extCount];
	m_extIndex = new double[m_extCount];
	}
else {
	m_extSrcArr = 0;
	m_extFlux = 0;
	m_extIndex = 0;
	}

cout << "Building " << m_extCount*m_mapCount << " extended source";
// zzz int count = m_extCount*m_mapCount;
if (m_extCount*m_mapCount>1)
	cout << "s";
cout << endl;
bool compatError = false;
for (int s=0; s<m_extCount; ++s) {
	double flux = m_extData.Flux(s);
	double index = m_extData.Index(s);
	FixFlag fixFlag = AllFixed;
	if (flux<0) {
		flux = -flux;
		fixFlag = fixFlag | FluxFree;
		}
	if (index<0) {
		index = -index;
		fixFlag = fixFlag | IndexFree;
		}
	for (int m=0; m<m_mapCount; ++m) {
		AlikeExtMap& extMap = m_extSrcArr[ExtIndex(m, s)];
		extMap.BuildFrom(m_extData.ExtMap(m, s), m_energyInf, m_energySup, index);
		extMap.SetFixflag(fixFlag);
		extMap.SetCoeff(flux);
		if (!extMap.CompatibleWith(m_expMaps[m])) {
			compatError = true;
			cerr << "ERROR: Extended source file " << extMap.GetFileName() << " is not compatible with " << m_expMaps[m].GetFileName() << endl;
			}
		}
	}
return !compatError;
}

bool RoiMulti::SetLogfile(const char* fileName)
{
if (m_logFile)
	m_logFile->close();
m_logFile = new ofstream(fileName);
if (!m_logFile->is_open()) {
	cerr << "Failed opening " << fileName << " for output" << endl;
	delete m_logFile;
	m_logFile = 0;
	return false;
	}
m_logFile->rdbuf()->pubsetbuf(0,0);
return true;
}

void RoiMulti::CloseLogFile()
{
delete m_logFile;
m_logFile = 0;
}


/**
void RoiMulti::Log(const char* text)
{
if (m_logFile)
	(*m_logFile) << text << endl;
}
*/


int RoiMulti::DoFit(
	const SourceDataArray& srcDataArr,
	double ranal,
	double ulcl,
	double loccl,
	int chatter,
	const char* sourceCheckTS,
	double minSourceTS)
{
if (!m_mapCount) {
	cerr << "No maps loaded: Fitting not possible" << endl;
	return 2;
	}
char option1[8];
char option2[8];
if (gROOT->GetVersionInt() < 53000) {
	strcpy(option1,"LLNM");
	strcpy(option2,"LLNEM");
} else {
	strcpy(option1,"LNM");
	strcpy(option2,"LNEM");
}	
if (chatter==0) {
	strcat(option1, "Q");
	strcat(option2, "Q");
	}
else if (chatter!=1) {
	strcat(option1, "V");
	strcat(option2, "V");
	}

/// m_countsHist.SetDirectory(0);

SetSources(srcDataArr, ranal, ulcl, loccl);
InitParams();

return DoFit(option1, option2, sourceCheckTS, minSourceTS);
}


static double LimitLatitude(double b)
{
if (b<-90.0)
	return -90.0;
else if (b>90.0)
	return 90.0;
return b;
}

static void EvalLimits(const SourceData& srcData, SourceLimits& srcLimits)
{
if (srcData.loclimit==0.0)
	srcLimits.lLim = srcLimits.bLim = false;
else {
	double cosb = cos(srcData.srcB*DEG2RAD);
	double longVar = 360.0;
	if (cosb!=0.0)
		longVar = srcData.loclimit/cosb;
	if (longVar>=180.0)
		srcLimits.lLim = false;
	else {
		srcLimits.lLim = true;
		srcLimits.lMin = srcData.srcL-longVar;
		srcLimits.lMax = srcData.srcL+longVar;
		}
	srcLimits.bLim = true;
	srcLimits.bMin = LimitLatitude(srcData.srcB-srcData.loclimit);
	srcLimits.bMax = LimitLatitude(srcData.srcB+srcData.loclimit);
	}
}




void RoiMulti::SetSources(const SourceDataArray& srcDataArr, double ranal, double ulcl, double loccl)
{
m_ranal = ranal;
m_ulcl = ulcl;
m_loccl = loccl;
m_inSrcDataArr = srcDataArr;
m_inSrcDataArr.SortByFlux();
m_outSrcDataArr = m_inSrcDataArr;
m_srcCount = m_inSrcDataArr.Count();
delete[] m_sources;
delete[] m_srcLimits;
delete[] m_fitInfo;
m_sources = new AlikeSourceMap[m_srcCount*m_mapCount];
m_srcLimits = new SourceLimits[m_srcCount];
m_fitInfo = new FitInfo[m_srcCount];

/// Correcting and printing the source input values
if (m_srcCount)
	cout << endl << "Sources [name sqrt(minTS) flux index position(l, b)]" << endl;
else
	cout << endl << "No point sources considered" << endl;
for (int i=0; i<m_srcCount; ++i) {
	int& fixFlag = m_inSrcDataArr[i].fixflag;
	if (fixFlag & (PosFree | IndexFree))	/// If any of Pos or Index are free => Flux will be free too
		fixFlag = fixFlag | FluxFree;

	if (fixFlag & PosFree)
		EvalLimits(m_inSrcDataArr[i], m_srcLimits[i]);
	else
		m_srcLimits[i].lLim = m_srcLimits[i].bLim = false;

	double& index = m_inSrcDataArr[i].index;
	if (index<0)	/// Adopting a positive value of index by convention
		index = -index;

	if (m_srcCount>9 && i<9)
		cout << "0";
	cout << i+1 << " " << m_inSrcDataArr[i].label;
	cout << " " << sqrt(m_inSrcDataArr[i].minTS);
	cout << " " << m_inSrcDataArr[i].flux;
	if (fixFlag & FluxFree)
		cout << " free ";
	else
		cout << " fixed ";
	cout << index;
	if (fixFlag & IndexFree)
		cout << " free ";
	else
		cout << " fixed ";
	cout << "(" << m_inSrcDataArr[i].srcL << ", " << m_inSrcDataArr[i].srcB << ")";
	if (fixFlag & PosFree) {
		if (!m_srcLimits[i].lLim && !m_srcLimits[i].bLim)
			cout << " free";
		else {
			if (m_srcLimits[i].lLim)
				cout << ", l free in [" << m_srcLimits[i].lMin << ".." << m_srcLimits[i].lMax << "]";
			else
				cout << ", l free";
			if (m_srcLimits[i].bLim)
				cout << ", b free in [" << m_srcLimits[i].bMin << ".." << m_srcLimits[i].bMax << "]";
			else
				cout << ", b free";
			}
		}
	else
		cout << " fixed";
	cout << endl;
	}

int src = 0;
for (int exp=0; exp<m_mapCount; ++exp) {
	const AgileMap& expMap = m_expMaps[exp];
	double theta = m_thetaArr[exp];
	for (int i=0; i<m_srcCount; ++i) {
		double srcL = m_inSrcDataArr[i].srcL;
		double srcB = m_inSrcDataArr[i].srcB;
		double index = m_inSrcDataArr[i].index;
		double flux = m_inSrcDataArr[i].flux;
		Double_t exposure = EvalExposure(srcL, srcB, expMap);
		m_sources[src].Set(&m_psfTab, expMap, m_energyInf, m_energySup, theta, srcL, srcB, index, flux, exposure);

	double cnt = 0;
	for (int r=0; r<m_sources[src].Rows(); ++r)
		for (int c=0; c<m_sources[src].Cols(); ++c)
			cnt += m_sources[src](r,c);
	cout << "Total SRC = " << cnt << endl;

		m_sources[src].SetFixflag(m_inSrcDataArr[i].fixflag);
		m_sources[src].SetLabel(m_inSrcDataArr[i].label);
		m_sources[src].SetMinTS(m_inSrcDataArr[i].minTS);
		m_sources[src].SetULCL(ulcl*ulcl);
		m_sources[src].SetLocCL(loccl);
		++src;
		}
	}
}



double RoiMulti::GetTotalExposure(int source) const
{
double exposure = 0;
if (source<m_srcCount && source>=0)
	for (int map=0; map<m_mapCount; ++map)
		exposure +=  m_sources[map*m_srcCount+source].GetExp()*m_sources[map*m_srcCount+source].GetNormFactor();
return exposure;
}


Double_t RoiMulti::GetSourceTS(const char* label) const
{
for (int i=0; i<m_srcCount; ++i)
	if (m_sources[i].GetLabel()==label)
		return m_sources[i].GetTS();
return 0;
}



void RoiMulti::Move(double lCenter, double bCenter)
{
m_lCenter = lCenter;
m_bCenter = bCenter;


	/// zzz temporary
	double totCts = 0;
	for (int r=0; r<m_ctsMaps[0].Rows(); ++r)
		for (int c=0; c<m_ctsMaps[0].Cols(); ++c)
			totCts += m_ctsMaps[0](r,c);
	cout << "Total CTS IN MOVE = " << totCts << endl;


AlikeCircle* masks = new AlikeCircle[m_alignedMaps?1:m_mapCount];
if (m_alignedMaps) {	/// Create a single mask if all the maps are aligned
	masks[0].Reset(lCenter, bCenter, m_ranal, m_ctsMaps[0]);
	m_binCount = masks[0].GetNoElements() * m_mapCount;
	}
else {	/// Create a mask per counts map otherwise
	m_binCount = 0;
	for (int cts=0; cts<m_mapCount; ++cts) {
		const AgileMap& map = m_ctsMaps[cts];
		masks[cts].Reset(lCenter, bCenter, m_ranal, map);
		m_binCount += masks[cts].GetNoElements();
		}
	}

delete[] m_binList;
delete[] m_binAddr;
m_binList = new int[m_binCount];
m_binAddr = new int[m_binCount];

m_countsHist.SetBins(m_binCount, -0.5, m_binCount-0.5);

int binProg = 0;
for (int cts=0; cts<m_mapCount; ++cts) {
	const AlikeCircle& mask = masks[m_alignedMaps?0:cts];
	int binCount = mask.GetNoElements();
	const AgileMap& map = m_ctsMaps[cts];
	for (int j=0; j<binCount; ++j) {
		m_binList[binProg] = mask.GetBinNo(j);
		m_binAddr[binProg] = cts;
		m_countsHist.SetBinContent(binProg+1, map(m_binList[binProg]));
		++binProg;
		}
	}
delete[] masks;
}


/**
static void PrintParams(int count, const double* params)
{
cout << "Params[" << count << "]:";
for (int i=0; i<count; ++i)
	cout << " " << params[i];
cout << endl;
}
*/

const char c_progChar[] = { '-', '\\', '|', '/' };
static int c_progIndx = 0;



Double_t RoiMulti::FitFunction(Double_t *x, Double_t *par)
{

/// m_fitParams.UpdateValues(par); /// zzz maybe this can be done when bin==0 only

int bin = int(x[0]);
if (!bin) {
	++s_fitterIterations;
	if (s_fitterIterations%100==0) {
		cout << c_progChar[c_progIndx] << '\r' << flush;
		c_progIndx = (c_progIndx+1)%4;
		}
	}
int map = m_binAddr[bin];
int srcOfs = map*m_srcCount;
const AlikeDiffMap& galSrc = m_galSrc[map];
const AlikeDiffMap& isoSrc = m_isoSrc[map];

double galCoeff = par[GalCoeffPar(map)];
if (m_galMode==DiffTogether)
	galCoeff *= m_galCoeffs[map];
double isoCoeff = par[IsoCoeffPar(map)];
if (m_isoMode==DiffTogether)
	isoCoeff *= m_isoCoeffs[map];

double result = galCoeff * galSrc(m_binList[bin]);
result += isoCoeff * isoSrc(m_binList[bin]); ///  * isoSrc.GetNormFactor();
for (int s=0; s<m_extCount; ++s)
	result += par[ExtCoeffPar(s)] * m_extSrcArr[ExtIndex(map,s)](m_binList[bin]);

for (int s=0; s<m_srcCount; ++s) {
	AlikeSourceMap& source = m_sources[srcOfs+s];
	AlikePsfSource::Changes changes = source.SetSrcData(par[SrcLPar(s)], par[SrcBPar(s)], par[SrcIdxPar(s)]);
	if (s_fitterChanges) {
		if (changes & AlikePsfSource::IndexChanged)
			s_fitterChanges[s*2] += 1;
		if (changes & AlikePsfSource::PositionChanged)
			s_fitterChanges[s*2+1] += 1;
		}
	result += par[SrcFluxPar(s)] * source(m_binList[bin]) * source.GetNormFactor() * source.GetExp() * m_fluxScaleFactorMulInv;
	}

/**
if (bin==0 && s_fitterIterations<0)
	PrintParams(m_sourceParOffset + 4*m_srcCount, par);
*/
return result;
}




void RoiMulti::GetExtPar(int source)
{
double param = m_model.GetParameter(ExtCoeffPar(source));
for (int map=0; map<m_mapCount; ++map)
	m_extSrcArr[ExtIndex(map, source)].SetCoeff(param);
}


void RoiMulti::Nullify(int source)
{
m_model.FixParameter(SrcFluxPar(source), 0);
m_model.FixParameter(SrcLPar(source), m_model.GetParameter(SrcLPar(source)));
m_model.FixParameter(SrcBPar(source), m_model.GetParameter(SrcBPar(source)));
m_model.FixParameter(SrcIdxPar(source), m_model.GetParameter(SrcIdxPar(source)));
} 


void RoiMulti::SetSrcPars(int source)
{
SetSrcFlux(source);
m_model.SetParameter(SrcLPar(source), m_sources[source].GetSrcL());
m_model.SetParameter(SrcBPar(source), m_sources[source].GetSrcB());
m_model.SetParameter(SrcIdxPar(source), m_sources[source].GetIndex());
}




void RoiMulti::ReleaseSrcPos(int source)
{
const SourceLimits& srcLimits = m_srcLimits[source];
if (srcLimits.lLim)
	m_model.SetParLimits(SrcLPar(source), srcLimits.lMin, srcLimits.lMax);
else
	m_model.ReleaseParameter(SrcLPar(source));
if (srcLimits.bLim)
	m_model.SetParLimits(SrcBPar(source), srcLimits.bMin, srcLimits.bMax);
else
	m_model.ReleaseParameter(SrcBPar(source));
}


void RoiMulti::FixSrcPos(int source, double l, double b)
{
m_model.FixParameter(SrcLPar(source), l);
m_model.FixParameter(SrcBPar(source), b);
}


void RoiMulti::SetSrcPos(int source, double l, double b)
{
m_model.SetParameter(SrcLPar(source), l);
m_model.SetParameter(SrcBPar(source), b);
}



void RoiMulti::InitParams()
{
/**
The fixFlag is here changed according to the parameters galMode and isoMode:
DiffFixed    = 0. All the coeffs are kept constant at the specified value
DiffDefault  = 1. The fix flags are not changed
DiffFree     = 2. All the coeffs are variable
DiffTogether = 3. All the coeffs are proportionally variable (single parameter)
*/



m_fluxScaleFactor = FLUX_SCALE_FACTOR;
m_fluxScaleFactorMulInv = 1.0 / FLUX_SCALE_FACTOR;
m_fluxUpperBound =  FLUX_UPPER_BOUND;

cout << endl << "Flux factor: " << m_fluxScaleFactor << "; upper bound: " << m_fluxUpperBound << endl;

m_galParCount = m_isoParCount = m_mapCount;
if (m_galMode==DiffTogether)
	m_galParCount = 1;
if (m_isoMode==DiffTogether)
	m_isoParCount = 1;
m_diffParCount = m_galParCount+m_isoParCount;
m_sourceParOffset = m_diffParCount+m_extCount;

m_diffPars[0].ResizeTo(m_diffParCount, m_srcCount, 2); /// Gal
m_diffPars[1].ResizeTo(m_diffParCount, m_srcCount, 2); /// Iso
m_diffErrs[0].ResizeTo(m_diffParCount, m_srcCount, 2); /// Gal
m_diffErrs[1].ResizeTo(m_diffParCount, m_srcCount, 2); /// Iso

m_model = TF1("Model", SingleFitFunction, 0, 0, m_sourceParOffset + 4*m_srcCount);

/// m_fitParams.Init(m_sourceParOffset + 4*m_srcCount);


/// In case diffMode==DiffTogether the diff parameter is set to 1
/// and the FitFunction will multiply by the original coefficient
cout << endl << "Parameters" << endl;
char parname[1024];
if (m_galMode==DiffTogether) {
	m_model.SetParName(GalCoeffPar(0), "Galactic");
	SetGalPar(0, 1);
	ReleaseGalPar(0);
	cout << GalCoeffPar(0)+1 << " - Galactic:";
	const char* sep[2] = { m_mapCount>1?" [":" ", ", " };
	for (int map=0; map<m_mapCount; ++map) {
		m_galSrc[map].SetFixflag(FluxFree);
		cout << sep[bool(map)] << m_galCoeffs[map];
		}
	if (m_mapCount>1)
		cout << "] free together" << endl;
	else
		cout << " free" << endl;
	}
else
	for (int map=0; map<m_mapCount; ++map) {
		sprintf(parname, "Galactic_%i", map+1);
		if (m_mapCount==1)
			strcpy(parname, "Galactic");
		m_model.SetParName(GalCoeffPar(map), parname);
		AlikeDiffMap& galSrc = m_galSrc[map];
		if (m_galMode==0)
			galSrc.SetFixflag(AllFixed);
		else if (m_galMode==2)
			galSrc.SetFixflag(FluxFree);
		double coeff = galSrc.GetCoeff();
		cout << GalCoeffPar(map)+1 << " - " << parname << ": " << coeff;
		if (galSrc.GetFixflag() & FluxFree) {
			SetGalPar(map);
			ReleaseGalPar(map);
			cout << " free" << endl;
			}
		else {
			FixGalPar(map);
			cout << " fixed" << endl;
			}
		for (int i=0; i<m_srcCount; ++i) {
			SetDPM(false, map, i, false, coeff);
			SetDPM(false, map, i, true, coeff);
			}
		}

if (m_isoMode==DiffTogether) {
	m_model.SetParName(IsoCoeffPar(0), "Isotropic");
	SetIsoPar(0, 1);
	ReleaseIsoPar(0);
	cout << IsoCoeffPar(0)+1 << " - Isotropic:";
	const char* sep[2] = { m_mapCount<2?" ":" [", ", " };

	for (int map=0; map<m_mapCount; ++map) {
		m_isoSrc[map].SetFixflag(FluxFree);
		cout << sep[bool(map)] << m_isoCoeffs[map];
		}
	if (m_mapCount>1)
		cout << "] free together" << endl;
	else
		cout << " free" << endl;
	}
else 
	for (int map=0; map<m_mapCount; ++map) {
		sprintf(parname, "Isotropic_%i", map+1);
		if (m_mapCount==1)
			strcpy(parname, "Isotropic");
		m_model.SetParName(IsoCoeffPar(map), parname);
		AlikeDiffMap& isoSrc = m_isoSrc[map];
		if (m_isoMode==0)
			isoSrc.SetFixflag(AllFixed);
		else if (m_isoMode>=2)
			isoSrc.SetFixflag(FluxFree);
		double coeff = isoSrc.GetCoeff();
		cout << IsoCoeffPar(map)+1 << " - " << parname << ": " << coeff;
		if (isoSrc.GetFixflag() & FluxFree) {
			SetIsoPar(map);
			ReleaseIsoPar(map);
			cout << " free" << endl;
			}
		else {
			FixIsoPar(map);
			cout << " fixed" << endl;
			}
		for (int i=0; i<m_srcCount; ++i) {
			SetDPM(true, map, i, false, coeff);
			SetDPM(true, map, i, true, coeff);
			}
		}


for (int i=0; i<m_extCount; ++i) {
	m_model.SetParName(ExtCoeffPar(i), m_extData.Name(i));
	cout << ExtCoeffPar(i) << " - " << m_extData.Name(i) << ": " << fabs(m_extData.Flux(i));
	if (m_extSrcArr[ExtIndex(0, i)].GetFixflag()) {
		SetExtPar(i);
		ReleaseExtPar(i);
		cout << " free" << endl;
		}
	else {
		FixExtPar(i);
		cout << " fixed" << endl;
		}
	}


for (int i=0; i<m_srcCount; ++i) {
	SetSrcPars(i);

	sprintf(parname, "%s_f", m_sources[i].GetLabel().c_str());
	m_model.SetParName(SrcFluxPar(i), parname);
	cout << SrcFluxPar(i)+1 << " - " << parname << ": " << m_sources[i].GetFlux();
	if (m_sources[i].GetFixflag() & FluxFree) {
		ReleaseSrcFlux(i);
		cout << " free" << endl;
		}
	else {
		FixSrcFlux(i);
		cout << " fixed" << endl;
		}

	sprintf(parname, "%s_x", m_sources[i].GetLabel().c_str());
	m_model.SetParName(SrcIdxPar(i), parname);
	cout << SrcIdxPar(i)+1 << " - " << parname << ": " << m_sources[i].GetIndex();
	if (m_sources[i].GetFixflag() & IndexFree) {
		ReleaseSrcIndex(i);
		cout << " free" << endl;
		}
	else {
		FixSrcIndex(i);
		cout << " fixed" << endl;
		}

	sprintf(parname, "%s_l", m_sources[i].GetLabel().c_str());
	m_model.SetParName(SrcLPar(i), parname);
	cout << SrcLPar(i)+1 << " - " << parname << ": " << m_sources[i].GetSrcL();
	if (m_sources[i].GetFixflag() & PosFree) {
		ReleaseSrcPos(i);
		cout << " free" << endl;
		}
	else {
		FixSrcPos(i);
		cout << " fixed" << endl;
		}

	sprintf(parname, "%s_b", m_sources[i].GetLabel().c_str());
	m_model.SetParName(SrcBPar(i), parname);
	cout << SrcBPar(i)+1 << " - " << parname << ": " << m_sources[i].GetSrcB();
	if (m_sources[i].GetFixflag() & PosFree)
		cout << " free" << endl;
	else
		cout << " fixed" << endl;
	}
cout << endl;
}




void RoiMulti::GetSrcPars(int source)
{
double flux = PeekSrcFluxPar(source);
double l = m_model.GetParameter(SrcLPar(source));
double b = m_model.GetParameter(SrcBPar(source));
double idx = m_model.GetParameter(SrcIdxPar(source));
for (int cts=0; cts<m_mapCount; ++cts) {
	m_sources[m_srcCount*cts+source].SetFlux(flux>0?flux:0);
	m_sources[m_srcCount*cts+source].SetSrcData(l, b, idx>0?idx:0);
	}
}




void RoiMulti::GetDiffPars()
{
for (int map=0; map<m_mapCount; ++map) {
	double coeff = m_model.GetParameter(GalCoeffPar(map));
	if (m_galMode==DiffTogether)
		coeff *= m_galCoeffs[map];
	m_galSrc[map].SetCoeff(coeff);
	coeff = m_model.GetParameter(IsoCoeffPar(map));
	if (m_isoMode==DiffTogether)
		coeff *= m_isoCoeffs[map];
	m_isoSrc[map].SetCoeff(coeff);
	}

/// Extended sources
for (int i=0; i<m_extCount; ++i) {
	double coeff = m_model.GetParameter(ExtCoeffPar(i));
	for (int map=0; map<m_mapCount; ++map)
		m_extSrcArr[ExtIndex(map, i)].SetCoeff(coeff);
	}
}


void RoiMulti::GetDiffErrs()
{
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
Double_t plus, minus, parab, globcc;

for (int i=0; i<m_mapCount; ++i) {
	if (m_galMode!=DiffTogether || i==0)
		fitter->GetErrors(GalCoeffPar(i), plus, minus, parab, globcc);
	double coeff = 1.0;
	if (m_galMode==DiffTogether)
		coeff = m_galCoeffs[i];
	AlikeDiffMap& galSrc = m_galSrc[i];
	galSrc.SetCoefferr(parab*coeff);
	galSrc.SetCoefflo(minus*coeff);
	galSrc.SetCoeffhi(plus*coeff);
	}

for (int i=0; i<m_mapCount; ++i) {
	if (m_isoMode!=DiffTogether || i==0)
		fitter->GetErrors(IsoCoeffPar(i), plus, minus, parab, globcc);
	double coeff = 1.0;
	if (m_isoMode==DiffTogether)
		coeff = m_isoCoeffs[i];
	AlikeDiffMap& isoSrc = m_isoSrc[i];
	isoSrc.SetCoefferr(parab*coeff);
	isoSrc.SetCoefflo(minus*coeff);
	isoSrc.SetCoeffhi(plus*coeff);
	}


/// Extended sources
for (int i=0; i<m_extCount; ++i) {
	fitter->GetErrors(ExtCoeffPar(i), plus, minus, parab, globcc);
	for (int map=0; map<m_mapCount; ++map) {
		AlikeDiffMap& extSrc = m_extSrcArr[ExtIndex(map, i)];
		extSrc.SetCoefferr(parab);
		extSrc.SetCoefflo(minus);
		extSrc.SetCoeffhi(plus);
		}
cerr << "Setting Ext source " << i+1 << " errors: parab=" << parab << ", minus=" << minus << ", plus=" << plus << endl;
	}
}


void RoiMulti::PrintDiffData()
{
cout << endl << "Gal sources [parNo coeff err lo hi]" << endl;
for (int i=0; i<m_mapCount; ++i) {
	AlikeDiffMap& galSrc = m_galSrc[i];
	cout << GalCoeffPar(i)+1 << " " << galSrc.GetCoeff() /// *m_galCoeff
		<< " " << galSrc.GetCoefferr() << " " << galSrc.GetCoefflo()
		<< " " << galSrc.GetCoeffhi() << endl;
	}
cout << "Iso sources [parNo coeff err lo hi]" << endl;
for (int i=0; i<m_mapCount; ++i) {
	AlikeDiffMap& isoSrc = m_isoSrc[i];
	cout << IsoCoeffPar(i)+1 << " " << isoSrc.GetCoeff() /// *m_isoCoeff
		<< " " << isoSrc.GetCoefferr() << " " << isoSrc.GetCoefflo()
		<< " " << isoSrc.GetCoeffhi() << endl;
	}
}

void RoiMulti::GetFluxErrors(int source)	/// zzz Add the errors about Index?
{
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
Double_t flux_plus, flux_minus, flux_parab, globcc;
fitter->GetErrors(SrcFluxPar(source), flux_plus, flux_minus, flux_parab, globcc);
for (int cts=0; cts<m_mapCount; ++cts) {
	m_sources[m_srcCount*cts+source].SetFluxerr(flux_parab*m_fluxScaleFactorMulInv);
	m_sources[m_srcCount*cts+source].SetFluxlo(flux_minus*m_fluxScaleFactorMulInv);
	m_sources[m_srcCount*cts+source].SetFluxhi(flux_plus*m_fluxScaleFactorMulInv);
	}
}


void RoiMulti::FitIndex(int source)
{
ReleaseSrcIndex(source);
SetSrcPars(source);
cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Finding best index" << endl;
if (gROOT->GetVersionInt() < 53000) {
	Fit("LLONEM", source, false, 1);
} else {
	Fit("LONEM", source, false, 1);
}
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
Double_t index_plus, index_minus, index_parab, globcc;
fitter->GetErrors(SrcIdxPar(source), index_plus, index_minus, index_parab, globcc);
for (int cts=0; cts<m_mapCount; ++cts)
	m_sources[m_srcCount*cts+source].SetIndexerr(index_parab);
GetSrcPars(source);
FixSrcIndex(source);
}




static Double_t s_oldErrorDef = 1;

static Double_t GetErrorDef()
{
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
return fitter->GetErrorDef();
}

static Double_t SetErrorDef(Double_t errDef)
{
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
s_oldErrorDef = fitter->GetErrorDef();
fitter->SetErrorDef(errDef);
return s_oldErrorDef;
}



void RoiMulti::GetSrcUL(int source)
{
ReleaseSrcFlux(source);
SetSrcPars(source);

Double_t newerrdef = m_sources[source].GetULCL();
if (gROOT->GetVersionInt() >= 53000) {newerrdef *= 0.5;}
Double_t olderrdef = SetErrorDef(newerrdef);
/// zzz
cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Finding upper level" << endl;
if (gROOT->GetVersionInt() < 53000) {
	Fit("LLONEM", source, false, 1);
} else {
	Fit("LONEM", source, false, 1);
}
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
Double_t flux = PeekSrcFluxPar(source);
Double_t flux_plus, flux_minus, flux_parab, globcc;
fitter->GetErrors(SrcFluxPar(source), flux_plus, flux_minus, flux_parab, globcc);
for (int cts=0; cts<m_mapCount; ++cts) {
	if (flux > 0.0)
		m_sources[m_srcCount*cts+source].SetFluxul(flux + flux_plus*m_fluxScaleFactorMulInv);
	else
		m_sources[m_srcCount*cts+source].SetFluxul(flux_plus*m_fluxScaleFactorMulInv - flux);
	}
SetErrorDef(olderrdef);
}





int GraphToEllipse(const TGraph* graph, AlikeSourceMap& source)
{
int sideCount = 0;
if (graph)
	sideCount = graph->GetN();
if (sideCount>=3) {
	cout << source.GetLabel() << ": Finding the best ellipse" << endl;
	Polygon polygon(sideCount);
	for (int i=0; i<sideCount; ++i)
		polygon[i] = Point(graph->GetX()[i], graph->GetY()[i]);
	double area;
	Point barycenter;
	polygon.GetBary(area, barycenter);
	double deform = cos(barycenter.y*DEG2RAD);
	Polygon defPoly(polygon);
	for (int i=0; i<sideCount; ++i) {
		Point& p = defPoly[i];
		p.x = barycenter.x+(p.x-barycenter.x)*deform;
		}
	area *= sqrt(deform);
	double radius = sqrt(area/M_PI);
	Ellipse ellipse;
	int result = PolyToEllipse(defPoly, barycenter, radius, ellipse);
	source.SetContour(polygon, ellipse, radius);
	return result;
	}
else {
	source.SetContour(Polygon(), Ellipse(), 0);
	cout << source.GetLabel() << ": Not enough points to find the best ellipse" << endl;
	}
return -1;
}



static void GetTsStat(Double_t& amin, Double_t& edm, int& iter)
{
Double_t errdef;
Int_t nvpar, nparx;
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
fitter->GetStats(amin, edm, errdef, nvpar, nparx);
iter = RoiMulti::GetFitterIterations();
}



static Int_t GetCovarStat()
{
Double_t amin, edm, errdef;
Int_t nvpar, nparx;
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
return fitter->GetStats(amin, edm, errdef, nvpar, nparx);
}


static double SqrtTS(double logL0, double logL1)
{
if (logL0>logL1)
	return sqrt(2.0*(logL0-logL1));
else
	return -sqrt(2.0*(logL1-logL0));
}

static double LnFact(double d)
{
double lf = 0;
if (d>1.0)
	for (double x=d; x>1; --x)
		lf += log(x);
return lf;
}


double RoiMulti::EvalFitFunction(double* params, double* data)
{
double like = 0;
if (data)
	*data = 0;
int zeroCount = 0;
s_fitterIterations = -2;
for (int i = 0; i<m_binCount; ++i) {
	int map = m_binAddr[i];
	const AgileMap& ctsMap = m_ctsMaps[map];
	int bin = m_binList[i];
	double mval = ctsMap(bin);
	double x = i;
	double fval = FitFunction(&x, params);
	if (fval)
		like -= mval*log(fval) - fval;
	else
		++zeroCount;
	if (data)
		*data += LnFact(mval);
	}
s_fitterIterations = 0;
if (zeroCount) {
	cerr << endl << "Error: EvalFitFunction returning zero " << zeroCount << " times over " << m_binCount << " bins" << endl;
	return 0;
	}
return like;
}


double RoiMulti::Likelihood(double* data)
{
int count = m_galParCount + m_isoParCount + 4*m_srcCount;
double* par = new double[count];
for (int p=0; p<count; ++p)
	par[p] = m_model.GetParameter(p);
double like = EvalFitFunction(par, data);
delete[] par;
return like;
}




void RoiMulti::MakeTSArrays(double*& tsVarArr, double*& tsConstArr)
{
tsVarArr = new double[m_srcCount];
tsConstArr = new double[m_srcCount];

int count = m_galParCount + m_isoParCount + 4*m_srcCount;
double* par = new double[count];

for (int s=0; s<m_srcCount; ++s) {
	const AlikeSourceMap& source = m_sources[s];
	if (source.GetFixflag()==AllFixed)
		tsConstArr[s] = tsVarArr[s] = 0;
	else {
		double circL = source.GetSrcL();
		double circB = source.GetSrcB();
		Move(circL, circB);
		/// Case with source
		for (int p=0; p<m_galParCount; ++p)
			par[p] = GetDPM(false, p, s, false);
		for (int p=0; p<m_isoParCount; ++p)
			par[p+m_galParCount] = GetDPM(true, p, s, false);
		for (int i=0; i<m_srcCount; ++i) {
			const AlikeSourceMap& source = m_sources[i];
			int index = m_galParCount+m_isoParCount+i*4;
			par[index] = source.GetFlux()*m_fluxScaleFactor;
			par[index+1] = source.GetIndex();
			par[index+2] = source.GetSrcL();
			par[index+3] = source.GetSrcB();
			}
		double data;
		double likeWith = EvalFitFunction(par, &data);

		/// Case without source, diffuse parameters unchanged (as if they were fixed)
		par[m_galParCount+m_isoParCount+s*4] = 0;
		double likeWithoutConst = EvalFitFunction(par, &data);

		/// Case without source, diffuse parameters variable
		for (int p=0; p<m_galParCount; ++p)
			par[p] = GetDPM(false, p, s, true);
		for (int p=0; p<m_isoParCount; ++p)
			par[p+m_galParCount] = GetDPM(true, p, s, true);
		double likeWithoutVar = EvalFitFunction(par, &data);

		tsConstArr[s] = SqrtTS(likeWithoutConst, likeWith);
		tsVarArr[s] = SqrtTS(likeWithoutVar, likeWith);
		}
	}
delete[] par;
}



/// This function is called during the first fitting loop if (Fixflag==FluxFree)
void RoiMulti::FitFlux(int source, const char* fitOpt)
{
ReleaseSrcFlux(source);
SetSrcPars(source);
cout << endl << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux free, position fixed" << endl;
Double_t amin1;
Fit(fitOpt, source, false, 1, &amin1);

/// double like1 = Likelihood();

for (int i=0; i<=source; ++i)
	GetSrcPars(i);
Nullify(source);
cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux=0, position fixed" << endl;
Double_t amin0;
Fit(fitOpt, source, true, 1, &amin0);

m_sources[source].SetTS(amin0 - amin1);
m_sources[source].PrintSqrtTS();

/**
double like0 = Likelihood();
cerr << endl << "amin1=" << amin1 << ", amin0=" << amin0 << ", L1=" << like1  << ", L0=" << like0 << endl;
cerr << "sqrt(TS)=" << sqrt(amin0 - amin1) << ", sqrt(TS2)=" << SqrtTS(like0, like1) << endl;
*/
if (m_sources[source].GetTS() < m_sources[source].GetMinTS() || (m_sources[source].GetFlux() < 0)) {
	Nullify(source);
	for (int i=0; i<=source; ++i)
		GetSrcPars(i);
	m_sources[source].SetFixflag(AllFixed);
	}
else {
	ReleaseSrcFlux(source);
	SetSrcPars(source);
	}
}


/// This function is called during the first fitting loop if (Fixflag&(PosFree|IndexFree))
void RoiMulti::FitPositionIndex(int source, const char* fitOpt)
{
ReleaseSrcFlux(source);
for (int i=0; i<SrcCount(); ++i) {
	if (m_sources[i].GetFixflag()) {
		FixSrcPos(i);
		FixSrcIndex(i);
		}
	}
SetSrcPars(source);
cout << endl << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux free, position and index fixed " << endl;
Double_t amin1;
Fit(fitOpt, source, false, 1, &amin1);

/// double like1 = Likelihood();

GetSrcPars(source);
Nullify(source);
cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux=0 " << endl;
Double_t amin0;
Fit(fitOpt, source, true, 1, &amin0);

m_sources[source].SetTS(amin0 - amin1);
m_sources[source].PrintSqrtTS();

/**
double like0 = Likelihood();
cerr << endl << "amin1=" << amin1 << ", amin0=" << amin0 << ", L1=" << like1  << ", L0=" << like0 << endl;
cerr << "sqrt(TS)=" << sqrt(amin0 - amin1) << ", sqrt(TS2)=" << SqrtTS(like0, like1) << endl;
*/

if (m_sources[source].GetTS()>m_sources[source].GetLocCL() && m_sources[source].GetFlux()>0) {
	cout << "Source " << source+1 << ": " << m_sources[source].GetLabel();
	SetSrcPars(source);
	ReleaseSrcFlux(source);
	int fitFlags = m_sources[source].GetFixflag();
	if ((fitFlags & PosFree) && (fitFlags & IndexFree)) {
		ReleaseSrcIndex(source);
		ReleaseSrcPos(source);
		cout << ": Flux, Index and Position free" << endl;
		}
	else if (fitFlags & PosFree) {
		ReleaseSrcPos(source);
		cout << ": Flux and Position free" << endl;
		}
	else if (fitFlags & IndexFree) {
		ReleaseSrcIndex(source);
		cout << ": Flux and Index free" << endl;
		}
	else
		cout << ": Flux free" << endl;
	double oldL = m_sources[source].GetSrcL();
	double oldB = m_sources[source].GetSrcB();
	double newL = oldL;
	double newB = oldB;
	int loopcount = 0;
	bool centergood = false;
	Double_t amin2 = amin1;
	do {
		cout << "Trying position (" << newL << ", " << newB << ")" << endl;
		SetSrcPos(source, newL, newB);
		Fit(fitOpt, source, false, fitFlags, &amin2);
		newL = m_model.GetParameter(SrcLPar(source));
		newB = m_model.GetParameter(SrcBPar(source));
		double dist = SphDistDeg(oldL, oldB, newL, newB);
		centergood = dist <= 0.25 * m_ranal;
		if (!centergood) {
			cout << "Analysis circle moved to (" << newL << "," << newB << ")" << endl;
			Move(newL, newB);
			oldL = newL;
			oldB = newB;
			}
		} while (++loopcount<1000 && !centergood);
	if (loopcount>1000)
		cout << "Location search failed to converge after " << loopcount << " iterations" << endl;
	Int_t result = GetCovarStat();
	cout << "Fit result = " << result << endl;
	if (amin2<amin1) {
		for (int i=0; i<=source; ++i)
			GetSrcPars(i);
		cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << " best position found in (l,b) = (";
		cout << m_sources[source].GetSrcL() << "," << m_sources[source].GetSrcB() << ")" << endl;
		}
	for (int i=0; i<=source; ++i) {
		FixSrcPos(i);
		FixSrcIndex(i);
		}
	GetDiffPars();

	cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux free, position fixed" << endl;
	Double_t amin1;
	Fit(fitOpt, source, false, 1, &amin1);

	/// double like1 = Likelihood();

	GetSrcPars(source);
	Nullify(source);
	cout << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux=0, position fixed" << endl;
	Double_t amin0;
	Fit(fitOpt, source, true, 0, &amin0);

	m_sources[source].SetTS(amin0 - amin1);
	m_sources[source].PrintSqrtTS();

	/**
	double like0 = Likelihood();
	cerr << endl << "amin1=" << amin1 << ", amin0=" << amin0 << ", L1=" << like1  << ", L0=" << like0 << endl;
	cerr << "sqrt(TS)=" << sqrt(amin0 - amin1) << ", sqrt(TS2)=" << SqrtTS(like0, like1) << endl;
	*/

	}
if (m_sources[source].GetTS() < m_sources[source].GetMinTS() || m_sources[source].GetFlux()<0) {
	Nullify(source);
	for (int i=0; i<=source; ++i)
		GetSrcPars(i);
	m_sources[source].SetFixflag(AllFixed);
	}
else {
	ReleaseSrcFlux(source);
	SetSrcFlux(source);
	FixSrcPos(source);
	FixSrcIndex(source);
	}
}



void RoiMulti::MakeCovarMat(int source)
{
TVirtualFitter* fitter = TVirtualFitter::GetFitter();
Int_t parCount = fitter->GetNumberTotalParameters();
TMatrixT<double> covar(parCount, parCount);
Double_t* covArr = fitter->GetCovarianceMatrix();
for (int i=0; i<parCount; ++i)
	for (int j=0; j<parCount; ++j)
		covar(i, j) = covArr[j+parCount*i];
TMatrixDEigen eigens(covar);
Double_t cosb = cosd(m_sources[source].GetSrcB());
// TMatrixD covar2(2, 2);
MatD covar2(2, 2, 2);
covar2(0,0) = sqrt(eigens.GetEigenValues()(SrcLPar(source), SrcLPar(source))*cosb*cosb);
double at2 = atan2d(
					eigens.GetEigenVectors()(SrcBPar(source), SrcLPar(source))*cosb,
					eigens.GetEigenVectors()(SrcLPar(source), SrcLPar(source))*cosb*cosb);

covar2(1,0) = at2;
covar2(0,1) = at2;
covar2(1,1) = sqrt(eigens.GetEigenValues()(SrcBPar(source), SrcBPar(source)));

/// m_sources[source].SetCovar(covar2);
/// cout << "Params " << SrcLPar(source)+1 << ", " << SrcBPar(source)+1 << endl;

/// zzz covar2.Print();
cout << "Covar Matrix" << endl;
cout << covar2(0,0) << " " << covar2(0,1) << endl;
cout << covar2(1,0) << " " << covar2(1,1) << endl;

}



void RoiMulti::MakeEllipse(int source)
{
cout << "Finding contour for source " << source+1 << endl;
s_fitterIterations = 0;
TGraph* graph = (TGraph*)gMinuit->Contour(40, SrcLPar(source), SrcBPar(source));
cout << "Done in " << s_fitterIterations-1 << " iterations" << endl;
GraphToEllipse(graph, m_sources[source]);
delete graph;
}


void RoiMulti::Loop2(const char* fitOpt)
{
cout << endl << "Begin second loop" << endl << endl;

bool diffParsStored = LoopExt(fitOpt);

double* savedFlux = new double[SrcCount()];
for (int source=0; source<SrcCount(); ++source) {
	savedFlux[source] = m_sources[source].GetFlux();
	Double_t fixedflux = m_sources[source].GetFlux();
	if (fixedflux == 0 || m_sources[source].GetFixflag()) {
		cout << endl<< "Considering source " << source+1 << ": " << m_sources[source].GetLabel() << endl;
		double srcL = m_sources[source].GetSrcL();
		double srcB = m_sources[source].GetSrcB();
		Move(srcL, srcB);
		for (int j=0; j<SrcCount(); ++j)
			if (!m_sources[j].GetFlux() || !m_sources[j].GetFixflag() || m_sources[j].DegDistance(srcL, srcB)>0.75*m_ranal)
				FixSrcFlux(j);
			else
				ReleaseSrcFlux(j);
		Nullify(source);
		cout << endl << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux=0, position fixed" << endl;
		Double_t amin0;
		Fit(fitOpt, source, true, 0, &amin0);
		if (!diffParsStored) {
			GetDiffPars();
			GetDiffErrs();
			/// PrintDiffData();
			}

		/// double like0 = Likelihood();


		ReleaseSrcFlux(source);
		SetSrcPars(source);
		cout << endl << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux free, position fixed" << endl;
		Double_t amin1;
		Fit(fitOpt, source, false, 1, &amin1);

		m_sources[source].SetTS(amin0 - amin1);
		m_sources[source].PrintSqrtTS();

		/**
		double like1 = Likelihood();
		cerr << endl << "amin0=" << amin0 << ", amin1=" << amin1 << ", L0=" << like0  << ", L1=" << like1 << endl;
		cerr << "sqrt(TS)=" << sqrt(amin0 - amin1) << ", sqrt(TS2)=" << SqrtTS(like0, like1) << endl;
		*/

		if (m_sources[source].GetFixflag() > FluxFree && m_sources[source].GetTS() > m_sources[source].GetLocCL()) {
			Double_t olderrdef = GetErrorDef();
			ReleaseSrcFlux(source);
			if (m_sources[source].GetFixflag() & PosFree) {
				ReleaseSrcPos(source);
				Double_t newerrdef = m_sources[source].GetLocCL();
				if (gROOT->GetVersionInt() >= 53000) {newerrdef *= 0.5;}
				SetErrorDef(newerrdef);
			}
			SetSrcPars(source);
			cout << endl << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux and position free, Errordef = " << m_sources[source].GetLocCL() << endl;
			Double_t amin2 = amin1;
			Fit(fitOpt, source, false, 3, &amin2);
			Int_t result = GetCovarStat();
			cout << "Fit result = " << result << endl;
			if (m_sources[source].GetFixflag() & PosFree) {
				if (amin2 < amin1) {
					GetSrcPars(source);
					srcL = m_sources[source].GetSrcL();
					srcB = m_sources[source].GetSrcB();
				}
//				MakeCovarMat(source);
				MakeEllipse(source);
				}
			FixSrcPos(source, srcL, srcB);
			SetErrorDef(olderrdef);
			cout << endl << "Source " << source+1 << ": " << m_sources[source].GetLabel() << ": Flux free, position fixed, Errordef = " << olderrdef << endl;
			Double_t amin1;
			Fit(fitOpt, source, false, 1, &amin1);

			m_sources[source].SetTS(amin0 - amin1);
			m_sources[source].PrintSqrtTS();

			/**
			double like1 = Likelihood();
			cerr << endl << "amin0=" << amin0 << ", amin1=" << amin1 << ", L0=" << like0  << ", L1=" << like1 << endl;
			cerr << "sqrt(TS)=" << sqrt(amin0 - amin1) << ", sqrt(TS2)=" << SqrtTS(like0, like1) << endl;
			*/

			}
		savedFlux[source] = PeekSrcFluxPar(source);
		GetSrcPars(source);

		if (m_sources[source].GetFixflag() & IndexFree)
			FitIndex(source);

		GetFluxErrors(source);
		if (!diffParsStored) {
			GetDiffPars();
			GetDiffErrs();
			/// PrintDiffData();
			diffParsStored = true;	/// Collecting the diff paramaters after the first (strongest) source
			}
		GetSrcUL(source);
		}
	else
		cout << "Skipping source " << source+1 << ": " << m_sources[source].GetLabel() << endl;

	if (fixedflux == 0) {
		Nullify(source);
		for (int cts=0; cts<m_mapCount; ++cts)
			m_sources[cts*m_srcCount+source].SetFlux(0);
		}
	else if (m_sources[source].GetFixflag() == 0)
		SetSrcFlux(source, fixedflux);
	else
		SetSrcPars(source);
	}

for (int source=0; source<SrcCount(); ++source)
	for (int cts=0; cts<m_mapCount; ++cts)
		m_sources[cts*m_srcCount+source].SetFlux(savedFlux[source]);

for (int source=0; source<SrcCount(); ++source) {
	double l = m_sources[source].GetSrcL();
	if (l>=360) {
		while (l>=360)
			l -= 360;
		double b = m_sources[source].GetSrcB();
		double x = m_sources[source].GetIndex();
		m_sources[source].SetSrcData(l, b, x);
		}
	}
delete[] savedFlux;


/** zzz Da rivedere con Andrew

for (int source=0; source<SrcCount(); ++source) {
	double index = m_sources[source].GetIndex();
	double srcL = m_sources[source].GetSrcL();
	double srcB = m_sources[source].GetSrcB();
	double exposure = 0;
	for (int cts=0; cts<m_mapCount; ++cts) {
		/// int mapIndex = cts*m_srcCount+source;
		double chFact = GetNormFactor(m_ctsMaps[cts].GetEmin(), m_ctsMaps[cts].GetEmax(), m_energyInf, m_energySup, index);
		double localExp = EvalExposure(srcL, srcB, m_expMaps[cts]);
		exposure += localExp*chFact;
		}
	/// Setting the exposure for the only source accessible from outside
	m_sources[source].SetExp(exposure);
	}
*/
}



static void SaveResults(SourceData& srcData, const AlikeSourceMap& map)
{
srcData.flux = map.GetFlux();
srcData.srcL = map.GetSrcL();
srcData.srcB = map.GetSrcB();
srcData.index = map.GetIndex();
srcData.TS = map.GetTS();
}



Int_t RoiMulti::Fit(const char* opt, int source, bool zero, int flags, Double_t* amin)
{
#ifdef DEBUG
    std::cout << "### pre-fitting parameters" << std::endl;
    m_model.Print("V");
#endif

/// zzz
for (int i=0; i<SrcCount(); ++i) {
	double cnt = 0;
	for (int r=0; r<m_sources[i].Rows(); ++r)
		for (int c=0; c<m_sources[i].Cols(); ++c)
			cnt += m_sources[i](r,c);
	cout << "Total CNT = " << cnt << ", NormFactor = " << m_sources[i].GetNormFactor() << ", GetExp() = " << m_sources[i].GetExp() << ", m_fluxScaleFactor = " << m_fluxScaleFactor << endl;
	}

if (source>=0)
	cout << "Fitting of source " << source+1;
else
	cout << "Fitting of extended source " << -source;
if (zero)
	cout << " (nullified)" << endl;
else {
	if (flags&1)
		cout << " Flux";
	if (flags&2)
		cout << " Position";
	if (flags&4)
		cout << " Index";
	if (flags&7)
		cout << " free";
	cout << endl;
	}

if (source>=0) { /// Exclude extended sources
	m_fitInfo[source].counts = 0;
	for (int i=0; i<m_binCount; ++i)
		m_fitInfo[source].counts += m_ctsMaps[m_binAddr[i]](m_binList[i]);
	}

s_fitterIterations = 0;
s_fitterChanges = new unsigned[m_srcCount*2];
for (int i=0; i<m_srcCount*2; ++i)
	s_fitterChanges[i] = 0;

Int_t fitResult = m_countsHist.Fit(&m_model, opt ,"", 0, m_binCount);

/// Evaluate the exposure again if the position was free


/// zzz
#define ADJUST_EXPOSURE
#ifdef ADJUST_EXPOSURE
if (source>=0) { /// Exclude extended sources
	if (flags&2) {
		for (int exp=0; exp<m_mapCount; ++exp) {
			double exposure = EvalExposure(m_sources[source].GetSrcL(), m_sources[source].GetSrcB(), m_expMaps[exp]);
			m_sources[exp*m_srcCount+source].SetExp(exposure);
			}
		}
	}
#endif

TVirtualFitter* fitter = TVirtualFitter::GetFitter();

/// Diffuse parameters and their errors are stored only when amin is requested too
/// zzz Check that this is correct
if (amin && source>=0) { /// zzz Extended sources not managed for now
	for (int par=0; par<m_galParCount; ++par)
		SetDPM(false, par, source, zero, m_model.GetParameter(GalCoeffPar(par)));
	for (int par=0; par<m_isoParCount; ++par)
		SetDPM(true, par, source, zero, m_model.GetParameter(IsoCoeffPar(par)));

	Double_t plus, minus, parab, globcc;
	for (int par=0; par<m_galParCount; ++par) {
		fitter->GetErrors(GalCoeffPar(par), plus, minus, parab, globcc);
		SetDPMErr(false, par, source, zero, parab);
		}
	for (int par=0; par<m_isoParCount; ++par) {
		fitter->GetErrors(IsoCoeffPar(par), plus, minus, parab, globcc);
		SetDPMErr(true, par, source, zero, globcc);
		}
	}


cout << "Fitting done in " << s_fitterIterations-1 << " iteration. Result=" << fitResult << endl;

int changesCount = 0;
int indexChanges = 0;
int positionChanges = 0;
for (int i=0; i<m_srcCount; ++i) {
	if (s_fitterChanges[i*2] || s_fitterChanges[i*2+1]) {
		++changesCount;
		cout << "Source " << i+1 << " changes: ";
		}
	if (s_fitterChanges[i*2]) {
		indexChanges += s_fitterChanges[i*2];
		cout << s_fitterChanges[i*2] << " times the index ";
		}
	if (s_fitterChanges[i*2+1]) {
		positionChanges += s_fitterChanges[i*2+1];
		cout << s_fitterChanges[i*2+1] << " times the position";
		}
	if (s_fitterChanges[i*2] || s_fitterChanges[i*2+1])
		cout << endl;
	}
/**
if (changesCount==0)
	cout << "No changes in index or position" << endl;
else if (changesCount>1) {
	cout << "Total changes in index: " << indexChanges << endl;
	cout << "Total changes in position: " << positionChanges << endl;
	}
*/
delete[] s_fitterChanges;
s_fitterChanges = 0;

Double_t fcn, edm;
int iter;
GetTsStat(fcn, edm, iter);

if (amin) {
	*amin = fcn;
	if (source>=0) { /// Exclude extended sources
		if (zero) {
			m_fitInfo[source].fcn0 = fcn;
			m_fitInfo[source].edm0 = edm;
			m_fitInfo[source].iter0 = iter;
			}
		else {
			m_fitInfo[source].fcn1 = fcn;
			m_fitInfo[source].edm1 = edm;
			m_fitInfo[source].iter1 = iter;
			}
		}
	}


if (m_logFile && source>=0) { /// Exclude extended sources for now
	/**
	double data;
	double likelihood = Likelihood(&data);
	*/
	
	ofstream& logFile = *m_logFile;
	char textFlags[4] = "ZFF";
	if (!zero) {
		textFlags[0] = flags&1 ? 'V' : 'F';
		textFlags[1] = flags&2 ? 'V' : 'F';
		textFlags[2] = flags&4 ? 'V' : 'F';
		textFlags[3] = 0;
		}
	if (source>=0)
		logFile << source+1 << " " << m_sources[source].GetLabel() << " " << textFlags;
	else
		logFile << -source << " " << m_extData.Name(-source-1) << " " << textFlags;

	char sep[2] = {' ', ','};
	int count = SingleGalPar() ? 1 : m_mapCount;
	for (int map=0; map<count; ++map)
		logFile << sep[bool(map)] << GetFinalDPM(false, map, source, zero);
	count = SingleIsoPar() ? 1 : m_mapCount;
	for (int map=0; map<count; ++map)
		logFile << sep[bool(map)] << GetFinalDPM(true, map, source, zero);	

	logFile << " " << m_lCenter << " " << m_bCenter << " " << fcn << " " << edm << " " << iter << endl;
	}

#ifdef DEBUG
    std::cout << "### post-fitting parameters" << std::endl;
    m_model.Print("V");
#endif

return fitResult;
}





void RoiMulti::FitExtendedSources(int source, const char* fitOpt, bool* diffParsStored)
{
ReleaseExtPar(source);
SetExtPar(source);
cout << endl << "Extended source " << source+1 << ": " << m_extData.Name(source) << ": Flux free" << endl;
Double_t amin0, amin1;
Fit(fitOpt, -source-1, false, 1, &amin1); /// zzz
if (!*diffParsStored) {
	*diffParsStored = true;
	GetDiffPars();
	GetDiffErrs();
	}

for (int i=0; i<=source; ++i)
	GetExtPar(i);
NullifyExt(source);
cout << "Extended source " << source+1 << ": " << m_extData.Name(source) << ": Flux=0" << endl;
Fit(fitOpt, -source-1, true, 1, &amin0); /// zzz

SetExtPar(source);
ReleaseExtPar(source);

for (int map=0; map<=m_mapCount; ++map)
	m_extSrcArr[ExtIndex(map, source)].SetTS(amin0 - amin1);
m_extSrcArr[ExtIndex(0, source)].PrintSqrtTS();
}




void RoiMulti::FixDistantSources(double L, double B, int source, bool extended)
{
int topSource = extended ? source : m_extCount;
for (int s=0; s<topSource; ++s)
	if (m_extSrcArr[s].GetFixflag()) {
		if (m_extSrcArr[s].DegDistance(L, B)>0.75*m_ranal)
			FixExtPar(s);
		else
			ReleaseExtPar(s);
		}
if (extended)
	return;
for (int s=0; s<source; ++s)
	if (m_sources[s].GetFixflag()) {
		if (m_sources[s].DegDistance(L, B)>0.75*m_ranal)
			FixSrcFlux(s);
		else
			ReleaseSrcFlux(s);
		}
}




bool RoiMulti::LoopExt(const char* fitOpt)
{
if (!m_extCount)
	return false;
bool diffParsStored = false;
for (int source=0; source<m_extCount; ++source) {
	if (m_extSrcArr[source].GetFixflag()) {
		cout << "Considering extended source " << source+1 << ": " << m_extData.Name(source) << endl;
		/**
		double circL = m_extSrcArr[source].GetMapCenterL();
		double circB = m_extSrcArr[source].GetMapCenterB();
		*/
		double circL = m_extSrcArr[source].GetBaryL();
		double circB = m_extSrcArr[source].GetBaryB();
		Move(circL, circB);
		FixDistantSources(circL, circB, source, true);
		FitExtendedSources(source, fitOpt, &diffParsStored);
		}
	else
		cout << "Skipping extended source " << source+1 << ": " << m_extData.Name(source) << endl;
	}
return diffParsStored;
}


void RoiMulti::Loop1(const char* fitOpt)
{
/// Move(m_sources[0].GetSrcL(), m_sources[0].GetSrcB());

/// zzz
for (int i=0; i<SrcCount(); ++i) {
	double cnt = 0;
	for (int r=0; r<m_sources[i].Rows(); ++r)
		for (int c=0; c<m_sources[i].Cols(); ++c)
			cnt += m_sources[i](r,c);
	cout << "Total CNT = " << cnt << endl;
	}


for (int i=0; i<m_extCount; ++i)
	if (m_extSrcArr[i].GetFixflag())
		NullifyExt(i);
for (int i=0; i<SrcCount(); ++i)
	if (m_sources[i].GetFixflag())
		Nullify(i);
cout << "Begin first loop" << endl << endl;


LoopExt(fitOpt);

for (int source=0; source<SrcCount(); ++source) {
	if (m_sources[source].GetFixflag()) {
		cout << "Considering source " << source+1 << ": " << m_sources[source].GetLabel() << endl;
		double circL = m_sources[source].GetSrcL();
		double circB = m_sources[source].GetSrcB();
		Move(circL, circB);
		FixDistantSources(circL, circB, source, false);
		if (m_sources[source].GetFixflag() & (PosFree|IndexFree))
			FitPositionIndex(source, fitOpt);
		else if (m_sources[source].GetFixflag() & FluxFree)
			FitFlux(source, fitOpt);
		}
	else
		cout << "Skipping source " << source+1 << ": " << m_sources[source].GetLabel() << endl;
	}
}



Int_t RoiMulti::DoFit(const char* option1, const char* option2, const char* sourceCheckTS, double minSourceTS)
{
Loop1(option1);
/**
Move(m_sources[0].GetSrcL(), m_sources[0].GetSrcB());

for (int i=0; i<SrcCount(); ++i)
	if (m_sources[i].GetFixflag())
		Nullify(i);
cout << "Begin first loop" << endl << endl;
for (int source=0; source<SrcCount(); ++source) {
	cout << "Considering source " << source+1 << ": " << m_sources[source].GetLabel() << endl;
	if (m_sources[source].GetFixflag()) {
		double circL = m_sources[source].GetSrcL();
		double circB = m_sources[source].GetSrcB();
		Move(circL, circB);
		for (int i=0; i<source; ++i)
			if (m_sources[i].DegDistance(circL, circB)>0.75*m_ranal)
				FixSrcFlux(i);
			else
				if (m_sources[i].GetFixflag())
					ReleaseSrcFlux(i);
		if (m_sources[source].GetFixflag() & (PosFree|IndexFree))
			FitPositionIndex(source, option1);
		else if (m_sources[source].GetFixflag() & FluxFree)
			FitFlux(source, option1);
		}
	}
*/

/// Skipping the second loop if the selected source has a low TS
if (sourceCheckTS) {
	bool skipSecondLoop = false;
	string sourceToCheck(sourceCheckTS);
	for (int s=0; s<SrcCount(); ++s)
		if (sourceToCheck==m_sources[s].GetLabel())
			if (m_sources[s].GetTS()<minSourceTS) {
				skipSecondLoop = true;
				break;
				}
	if (skipSecondLoop) {
		cout << "Source " << sourceCheckTS << " has a low TS, skipping the second loop" << endl;
		for (int s=0; s<SrcCount(); ++s) {
			string label = m_sources[s].GetLabel();
			/// m_outSrcDataArr[label].GetResultsFrom(m_sources[s]);
			SaveResults(m_outSrcDataArr[label], m_sources[s]);
			}
		return -1;
		}
	}

Loop2(option2);

for (int s=0; s<SrcCount(); ++s) {
	string label = m_sources[s].GetLabel();
	/// m_outSrcDataArr[label].GetResultsFrom(m_sources[s]);
	SaveResults(m_outSrcDataArr[label], m_sources[s]);
	}
cout << endl << "End of fitting." << endl;


PrintDiffParams();

return 0;
}


void RoiMulti::PrintDiffParams() const
{
bool diffVar = false; /// At least one diffuse source is variable
for (int map=0; map<m_mapCount && !diffVar; ++map)
	diffVar = m_galSrc[map].GetFixflag()!=AllFixed || m_isoSrc[map].GetFixflag()!=AllFixed;


if (diffVar && SrcCount()) {
	bool* varSrcArr = new bool[SrcCount()];
	bool varFlag = false;
	for (int s=0; s<SrcCount(); ++s) {
		varSrcArr[s] = m_inSrcDataArr[s].fixflag!=AllFixed;
		varFlag = varFlag || varSrcArr[s];
		}
	if (varFlag) {
		cout << endl << "Diffuse parameters tables" << endl;
		cout << "Values around each source, with(1) or without(0) it" << endl;
		size_t nameSize = strlen("Source");
		for (int s=0; s<SrcCount(); ++s)
			if (varSrcArr[s])
				if (GetSource(s).GetLabel().length()>nameSize)
					nameSize = GetSource(s).GetLabel().length();
		cout << "Coefficients" << endl;
		cout << setw(1+nameSize) << left << "Source" << right;
		for (int par=0; par<m_galParCount; ++par)
			cout << setw(10) << "gal1" << setw(10) << "gal0";
		for (int par=0; par<m_isoParCount; ++par)
			cout << setw(10) << "iso1" << setw(10) << "iso0";
		cout << endl;
		for (int s=0; s<SrcCount(); ++s)
			if (varSrcArr[s]) {
				cout << setw(1+nameSize) << left << GetSource(s).GetLabel() << right;
				for (int par=0; par<m_galParCount; ++par) {
					cout << setw(10) << GetFinalDPM(false, par, s, false);
					cout << setw(10) << GetFinalDPM(false, par, s, true);
					}
				for (int par=0; par<m_isoParCount; ++par) {
					cout << setw(10) << GetFinalDPM(true, par, s, false);
					cout << setw(10) << GetFinalDPM(true, par, s, true);
					}
				cout << endl;
				}
		cout << "Errors" << endl;
		cout << setw(1+nameSize) << left << "Source" << right;
		for (int par=0; par<m_galParCount; ++par)
			cout << setw(10) << "gal1" << setw(10) << "gal0";
		for (int par=0; par<m_isoParCount; ++par)
			cout << setw(10) << "iso1" << setw(10) << "iso0";
		cout << endl;
		for (int s=0; s<SrcCount(); ++s)
			if (varSrcArr[s]) {
				cout << setw(1+nameSize) << left << GetSource(s).GetLabel() << right;
				for (int par=0; par<m_galParCount; ++par) {
					cout << setw(10) << GetFinalDPMErr(false, par, s, false);
					cout << setw(10) << GetFinalDPMErr(false, par, s, true);
					}
				for (int par=0; par<m_isoParCount; ++par) {
					cout << setw(10) << GetFinalDPMErr(true, par, s, false);
					cout << setw(10) << GetFinalDPMErr(true, par, s, true);
					}
				cout << endl;
				}
		}
        delete []varSrcArr;
	}
}





double RoiMulti::GetFinalDPM(bool iso, int par, int source, bool zero) const
{
if (iso==true && m_isoMode==DiffTogether)
	return GetDPM(iso, 0, source, zero)*m_isoCoeffs[par];
else if (iso==false && m_galMode==DiffTogether)
	return GetDPM(iso, 0, source, zero)*m_galCoeffs[par];
return GetDPM(iso, par, source, zero);
}

double RoiMulti::GetFinalDPMErr(bool iso, int par, int source, bool zero) const
{
if (iso==true && m_isoMode==DiffTogether)
	return GetDPMErr(iso, 0, source, zero)*m_isoCoeffs[par];
else if (iso==false && m_galMode==DiffTogether)
	return GetDPMErr(iso, 0, source, zero)*m_galCoeffs[par];
return GetDPMErr(iso, par, source, zero);
}


bool RoiMulti::SingleGalPar() const
{
bool singleGal = m_galMode==DiffTogether || m_mapCount==1;
for (int map=1; map<m_mapCount && singleGal; ++map)
	if (m_galCoeffs[map]!=m_galCoeffs[0])
		singleGal = false;
return singleGal;
}

bool RoiMulti::SingleIsoPar() const
{
bool singleIso = m_isoMode==DiffTogether || m_mapCount==1;
for (int map=1; map<m_mapCount && singleIso; ++map)
	if (m_isoCoeffs[map]!=m_isoCoeffs[0])
		singleIso = false;
return singleIso;
}


void RoiMulti::Write(const char* fileName, bool skipFitInfo) const
{
if (!m_mapCount || !(m_srcCount+m_extCount)) {
	cerr << "Not writing file " << fileName << ": No data to write" << endl;
	return;
	}
//string srcoutname(string(fileName) + ".sources");
//ofstream output(srcoutname.c_str());
	ofstream output(fileName);

char parname[64];
output << "! DiffName, Flux, Err, +Err, -Err" << endl;

bool singleGal = SingleGalPar();
bool singleIso = SingleIsoPar();

if (singleGal) {
	AlikeDiffMap& gasMap = m_galSrc[0];
	output << "Galactic " << gasMap.GetCoeff() << " " << gasMap.GetCoefferr()
		<< " " << gasMap.GetCoeffhi() << " " << gasMap.GetCoefflo() << endl;
	}
else
	for (int map=0; map<m_mapCount; ++map) {
		sprintf(parname, "Galactic_%i", map+1);
		AlikeDiffMap& gasMap = m_galSrc[map];
		output << parname << " " << gasMap.GetCoeff() << " " << gasMap.GetCoefferr()
			<< " " << gasMap.GetCoeffhi() << " " << gasMap.GetCoefflo() << endl;
		}

if (singleIso) {
	AlikeDiffMap& isoMap = m_isoSrc[0];
	output << "Isotropic " << isoMap.GetCoeff() << " " << isoMap.GetCoefferr()
		<< " " << isoMap.GetCoeffhi() << " " << isoMap.GetCoefflo() << endl;
	}
else
	for (int map=0; map<m_mapCount; ++map) {
		sprintf(parname, "Isotropic_%i", map+1);
		AlikeDiffMap& isoMap = m_isoSrc[map];
		output << parname << " " << isoMap.GetCoeff() << " " << isoMap.GetCoefferr()
			<< " " << isoMap.GetCoeffhi() << " " << isoMap.GetCoefflo() << endl;
		}


if (ExtCount()) {
	output << "! ExtSrcName, sqrt(TS), Flux, Err, +Err, -Err" << endl;
	for (int i=0; i<ExtCount(); ++i) {
		const AlikeExtMap& extMap = m_extSrcArr[ExtIndex(0, i)];
		double BaryExposure = EvalExposure(extMap.GetBaryL(), extMap.GetBaryB(), m_expMaps[0]);
		output << m_extData.Name(i) << " " << sqrt(extMap.GetTS());
		output << " " << extMap.GetCoeff()/BaryExposure;
		output << " " << extMap.GetCoefferr()/BaryExposure;
		output << " " << extMap.GetCoeffhi()/BaryExposure;
		output << " " << extMap.GetCoefflo()/BaryExposure << endl;
		}
	}



if (m_srcCount) {
	output << "! SrcName, sqrt(TS), L_peak, B_peak, Counts, Err, Flux, Err, Index, Err";
	if (!skipFitInfo)
		output << " [fcn0 fcn1 edm0 edm1 iter0 iter1]";
	output << endl;
	
	for (int i=0; i<m_srcCount; ++i) {
	
		double exposure = GetTotalExposure(i);
	
		output << m_sources[i].GetLabel()
					<< " " << sqrt(m_sources[i].GetTS())
					<< " " << m_sources[i].GetSrcL()
					<< " " << m_sources[i].GetSrcB()
	/**
					<< " " << m_sources[i].GetCounts()
					<< " " << m_sources[i].GetCountsErr()
	*/
					<< " " << exposure*m_sources[i].GetFlux()
					<< " " << exposure*m_sources[i].GetFluxerr()
	
					<< " " << m_sources[i].GetFlux()
					<< " " << m_sources[i].GetFluxerr()
					<< " " << m_sources[i].GetIndex()
					<< " " << m_sources[i].GetIndexerr();
		if (m_sources[i].GetFixflag() && !skipFitInfo)
			output << " " << m_fitInfo[i].fcn0
					<< " " << m_fitInfo[i].fcn1
					<< " " << m_fitInfo[i].edm0
					<< " " << m_fitInfo[i].edm1
					<< " " << m_fitInfo[i].iter0
					<< " " << m_fitInfo[i].iter1;
		output << endl;
		}
	}
output.close();
}

static bool OpenConditionalBin(ofstream& f, bool sameValues, int row, int count)
{
	if (!sameValues) {
		f << "<td>";
		return true;
	}
	else
		if (row==0) {
			f << "<td rowspan=" << count << ">";
			return true;
		}
	return false;
}


static void WriteIsoDate(ofstream& f, const char* isoDdate)
{
	char date[64];
	strncpy(date, isoDdate, 32);
	date[33]='\0';
	if (strlen(date)==19 && date[10]=='T') {
		date[10] = 0;
		f << date << "&nbsp;" << date+11;
	}
	else if (date[0]==0)
		f << "&nbsp;";
	else
		f << isoDdate;
}

static void WriteIsoDate2(ofstream& f, const char* isoDdate)
{
	char date[64];
	strncpy(date, isoDdate, 32);
	date[33]='\0';
	if (strlen(date)==19 && date[10]=='T') {
		date[10] = 0;
		f << date << "T" << date+11;
	}
	else if (date[0]==0)
		f << "T";
	else
		f << isoDdate;
}

void RoiMulti::WriteSources(const char* fileName, bool skipFixed, bool skipEllipses) const
{
for (int i=0; i<m_srcCount; ++i) {
	if (skipFixed && !m_inSrcDataArr[i].fixflag)
		continue;
	string srcoutname(string(fileName) + "_" + m_sources[i].GetLabel() + ".source");
	ofstream srcout(srcoutname.c_str());
	srcout << "! Label, Fix, index, UL conf. level, srcloc conf. level, start l, start b, start flux, [ lmin , lmax ], [ bmin, bmax ]" << endl;
	srcout << "! sqrt(TS)" << endl ;
	srcout << "! L_peak, B_peak, Dist from initial position" << endl;
	const Ellipse& ellipse = m_sources[i].GetEllipse();
	srcout << "! L, B, Dist from initial position, r, a, b, phi" << endl;
	srcout << "! Counts, Err, +Err, -Err, UL" << endl;
	srcout << "! Flux, Err, +Err, -Err, UL, Exp" << endl;
	srcout << "! Index, Err" << endl;
	srcout << "! cts, fcn0, fcn1, edm0, edm1, iter0, iter1" << endl;
	
	//if (m_inSrcDataArr[i].fixflag) {
		srcout << "! Gal coeffs and errs" << endl;
		srcout << "! Gal zero coeffs and errs" << endl;
		srcout << "! Iso coeffs and errs" << endl;
		srcout << "! Iso zero coeffs and errs" << endl;
	//	}
	srcout << "! Start date, end date" << endl;
	srcout << "! Emin..emax, fovmin..fovmax, albedo, binsize, expstep, phasecode" << endl;
	srcout << m_sources[i].GetLabel()
				<< " " << m_inSrcDataArr[i].fixflag	/// Original flag
				<< " " << m_inSrcDataArr[i].index	/// Original index
				<< " " << sqrt(m_sources[i].GetULCL())
				<< " " << m_sources[i].GetLocCL()
				<< " " << m_inSrcDataArr[i].srcL
				<< " " << m_inSrcDataArr[i].srcB
	<< " " << m_inSrcDataArr[i].flux;
	
	if (m_srcLimits[i].lLim)
		srcout <<  " [ " << m_srcLimits[i].lMin << " , " << m_srcLimits[i].lMax << " ] ";
	else
		srcout << " [ -1 , -1 ] ";
	
	if (m_srcLimits[i].bLim)
		srcout << " [ " << m_srcLimits[i].bMin << " , " << m_srcLimits[i].bMax << " ] ";
	else
		srcout << " [ -1 , -1 ] ";
	
		srcout	<< endl;
	/// zzz m_sources[i].PrintAbphi();
	double dist = SphDistDeg(m_sources[i].GetSrcL(), m_sources[i].GetSrcB(), m_inSrcDataArr[i].srcL, m_inSrcDataArr[i].srcB);
	
	srcout << sqrt(m_sources[i].GetTS()) << endl << m_sources[i].GetSrcL() << " " << m_sources[i].GetSrcB() << " " << dist << endl;

	if (ellipse.horAxis!=0 && ellipse.verAxis!=0) {
		double dist2 = SphDistDeg(ellipse.center.x, ellipse.center.y, m_inSrcDataArr[i].srcL, m_inSrcDataArr[i].srcB);
		srcout << ellipse.center.x << " " << ellipse.center.y << " " << dist2 << " " << m_sources[i].GetRadius() << " " << ellipse.horAxis << " " << ellipse.verAxis << " " << ellipse.attitude*RAD2DEG << " " << endl;
	} else
		srcout << "-1 -1 -1 -1 -1 -1 -1 " << endl;

	double exposure = GetTotalExposure(i);
/**
	srcout << m_sources[i].GetCounts()
				<< " " << m_sources[i].GetCountsErr()
				<< " " << m_sources[i].GetCountshi()
				<< " " << m_sources[i].GetCountslo()
				<< " " << m_sources[i].GetCountsul()
				<< endl;
*/
	srcout << exposure*m_sources[i].GetFlux()
				<< " " << exposure*m_sources[i].GetFluxerr()
				<< " " << exposure*m_sources[i].GetFluxhi()
				<< " " << exposure*m_sources[i].GetFluxlo()
				<< " " << exposure*m_sources[i].GetFluxul()
				<< endl;

	srcout << m_sources[i].GetFlux()
				<< " " << m_sources[i].GetFluxerr()
				<< " " << m_sources[i].GetFluxhi()
				<< " " << m_sources[i].GetFluxlo()
				<< " " << m_sources[i].GetFluxul()
				<< " " << exposure /// m_sources[i].GetExp()
				<< endl;
	srcout << m_sources[i].GetIndex() << " " << m_sources[i].GetIndexerr() << endl;
	
	//AB
	
	srcout << m_fitInfo[i].counts;
	
	srcout << " " << m_fitInfo[i].fcn0;
	srcout << " " << m_fitInfo[i].fcn1;
	srcout << " " << m_fitInfo[i].edm0;
	srcout << " " << m_fitInfo[i].edm1;
	srcout << " " << m_fitInfo[i].iter0;
	srcout << " " << m_fitInfo[i].iter1 << endl;

	if (m_inSrcDataArr[i].fixflag) {
		const char* sep[2] = {"", ","};
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPM(false, m, i, false);
		srcout << ' ';
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPMErr(false, m, i, false);
		srcout << endl;
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPM(false, m, i, true);
		srcout << ' ';
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPMErr(false, m, i, true);
		srcout << endl;

		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPM(true, m, i, false);
		srcout << ' ';
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPMErr(true, m, i, false);
		srcout << endl;
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPM(true, m, i, true);
		srcout << ' ';
		for (int m=0; m<m_mapCount; ++m)
			srcout << sep[bool(m)] << GetFinalDPMErr(true, m, i, true);
		srcout << endl;
	} else  {
		srcout << "0 0\n0 0\n0 0\n0 0\n";
	}
	bool sameDate = true;
	for (int map=0; map<m_mapCount; ++map) {
		const AgileMap& m = m_ctsMaps[map];
		WriteIsoDate2(srcout, m.GetStartDate());
		srcout << " ";
		WriteIsoDate2(srcout, m.GetEndDate());
		srcout << endl;
		break;
	}
	for (int map=0; map<m_mapCount; ++map) {
		const AgileMap& m = m_ctsMaps[map];
		srcout << m.GetEmin() << ".." << m.GetEmax();
		if(m_mapCount > 1 && map != m_mapCount - 1)
			srcout << ",";
	}
	srcout << " ";
	for (int map=0; map<m_mapCount; ++map) {
		const AgileMap& m = m_ctsMaps[map];
		srcout << m.GetFovMin() << ".." << m.GetFovMax();
		if(m_mapCount > 1 && map != m_mapCount - 1)
			srcout << ",";
	}
	srcout << " ";
	const AgileMap& m = m_ctsMaps[0];
	srcout << m.GetAlbedo() << " ";
	srcout << m.GetYbin() << " ";
	srcout << m_expMaps[0].GetStep() << " ";
	srcout << m.GetPhaseCode() << endl;
	/*
	if (OpenConditionalBin(htmlout, sameEnergy, map, m_mapCount))
		htmlout << m.GetEmin() << ".." << m.GetEmax() << "</td>";
	if (OpenConditionalBin(htmlout, sameFOV, map, m_mapCount))
		htmlout << m.GetFovMin() << ".." << m.GetFovMax() << "</td>";
	
	if (OpenConditionalBin(htmlout, sameCenter, map, m_mapCount))
		htmlout <<m.GetMapCenterL() << ",&nbsp;" << m.GetMapCenterB() << "</td>";
	WriteConditionalBin(htmlout, m.GetAlbedo(), sameAlbedo, map, m_mapCount);
	WriteConditionalBin(htmlout, m.GetYbin(), sameBinSize, map, m_mapCount);
	WriteConditionalBin(htmlout, m_expMaps[map].GetStep(), sameStep, map, m_mapCount);
	WriteConditionalBin(htmlout, m.GetPhaseCode(), samePhCode, map, m_mapCount);
	*/
	srcout.close();
	if (skipEllipses)
		continue;
	m_sources[i].WritePolygon(srcoutname);
	m_sources[i].WriteEllipse(srcoutname);
	}
}


void RoiMulti::LogSources(const char* fileName, int iterNum) const
{
for (int i=0; i<m_srcCount; ++i) {
	if (!m_inSrcDataArr[i].fixflag)
		continue;


/// NUMIT L B TS FLUX FLUXERR FLUXUL SPECTRAL_INDEX FIXFLAG MINTS R EXP CTS CTSERROR CTSUL TOTEXP TOTNCOUNTS FCN0 FCN1 EDM0 EDM1 ITER0 ITER1 GAL ISO



	string srcoutname(string(fileName) + "_" + m_sources[i].GetLabel());
	ofstream srcout(srcoutname.c_str(), ios::app);

	srcout << iterNum;
	srcout << " " << m_sources[i].GetSrcL();
	srcout << " " << m_sources[i].GetSrcB();
	srcout << " " << m_sources[i].GetTS();
	srcout << " " << m_sources[i].GetFlux();
	srcout << " " << m_sources[i].GetFluxerr();
	srcout << " " << m_sources[i].GetFluxul();
	srcout << " " << m_sources[i].GetIndex();
	srcout << " " << m_inSrcDataArr[i].fixflag; /// zzz why not from the sources array?
	srcout << " " << m_sources[i].GetMinTS();


	if (m_sources[i].GetRadius())
		srcout << " " << m_sources[i].GetRadius();
	else
		srcout << " -1";

	srcout << " " << m_sources[i].GetExp();

	double exposure = GetTotalExposure(i);

	srcout << " " << exposure*m_sources[i].GetFlux();
	srcout << " " << exposure*m_sources[i].GetFluxerr();
	srcout << " " << exposure*m_sources[i].GetFluxul();
	srcout << " " << exposure;
	srcout << " " << m_fitInfo[i].counts;

	srcout << " " << m_fitInfo[i].fcn0;
	srcout << " " << m_fitInfo[i].fcn1;
	srcout << " " << m_fitInfo[i].edm0;
	srcout << " " << m_fitInfo[i].edm1;
	srcout << " " << m_fitInfo[i].iter0;
	srcout << " " << m_fitInfo[i].iter1;
	
	char sep[2] = {' ', ','};
	int count = SingleGalPar() ? 1 : m_mapCount;
	for (int map=0; map<count; ++map)
		srcout << sep[bool(map)] << GetFinalDPM(false, map, i, false);
	count = SingleIsoPar() ? 1 : m_mapCount;
	for (int map=0; map<count; ++map)
		srcout << sep[bool(map)] << GetFinalDPM(true, map, i, false);

	srcout << endl;

	}
}


static void StripExtension(char* str, const char* ext1, const char* ext2)
{
char* p = strstr(str, ext1);
if (!p)
	p = strstr(str, ext2);
if (p)
	*p = 0;
}





static void WriteConditionalBin(ofstream& f, double value, bool sameValues, int row, int count)
{
if (OpenConditionalBin(f, sameValues, row, count))
	f << value << "</td>";
}


void RoiMulti::WriteHtml(const char* fileName, const char* suffix) const
{
if (!m_mapCount || !(m_srcCount+m_extCount)) {
	cerr << "Not writing file " << fileName << ": No data to write" << endl;
	return;
	}

string htmloutname(fileName);
htmloutname += suffix;
ofstream htmlout(htmloutname.c_str());

bool singleGal = SingleGalPar();
bool singleIso = SingleIsoPar();

htmlout << "<html>" << endl << "<head><title>" << fileName << "</title></head>" << endl << "<body>" << endl;
htmlout << "<h1>" << fileName << "</h1>" << endl;

time_t theTime;
time(&theTime);
htmlout << ROIMULTI_VERSION << " - " << ctime(&theTime) << endl;

htmlout << "<h2>Input</h2>" << endl;


htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
htmlout << "<tr><td>Psf</td><td>" << m_psfTab.PsfFileName() << "</td></tr>" << endl;
htmlout << "<tr><td>Raeff</td><td>" << m_psfTab.RaeffFileName() << "</td></tr>" << endl;
htmlout << "<tr><td>Edp</td><td>" << m_psfTab.EdpFileName() << "</td></tr>" << endl;
htmlout << "</table>" << endl << "<br>" << endl;

htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;

htmlout << "<tr><td>Gal Mode</td><td>Iso Mode</td><td>Radius</td><td>ulcl</td><td>loccl</td></tr>" << endl;
htmlout << "<tr><td>" << m_isoMode << "</td><td>" << m_galMode << "</td><td>" << m_ranal
			<< "</td><td>" << m_ulcl << "</td><td>" << m_loccl << "</td></tr>" << endl;
htmlout << "</table>" << endl << "<br>" << endl;


htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
htmlout << "<tr><td>Map</td><td>Name</td><td>Theta</td><td>GalCoeff</td><td>IsoCoeff</td></tr>" << endl;
for (int map=0; map<m_mapCount; ++map) {
	const AgileMap& m = m_ctsMaps[map];
	char mapName[256];
	strcpy(mapName, m.GetFileName());
	StripExtension(mapName, ".cts.gz", ".cts");
	htmlout << "<tr><td>" << map+1
		<< "</td><td>" << mapName
		<< "</td><td>" << m_thetaArr[map]
		<< "</td><td>" << m_galCoeffs[map] << "&nbsp;" << ((m_galSrc[map].GetFixflag() & FluxFree) ? "free":"fixed")
		<< "</td><td>" << m_isoCoeffs[map] << "&nbsp;" << ((m_isoSrc[map].GetFixflag() & FluxFree) ? "free":"fixed")
		<< "</td></tr>" << endl;
	}
htmlout << "</table>" << endl;
htmlout << "<br>" << endl;

bool sameDate = true;
bool sameEnergy = true;
bool sameFOV = true;
bool sameCenter = true;
bool sameAlbedo = true;
bool sameBinSize = true;
bool sameStep = true;
bool samePhCode = true;
for (int map=0; map<m_mapCount-1; ++map) {
	const AgileMap& tMap = m_ctsMaps[map];
	const AgileMap& nMap = m_ctsMaps[map+1];
	if (sameDate)
		if (strcmp(tMap.GetStartDate(), nMap.GetStartDate()) || strcmp(tMap.GetEndDate(), nMap.GetEndDate()))
			sameDate = false;
	if (sameEnergy)
		if (tMap.GetEmin()!=nMap.GetEmin() || tMap.GetEmax()!=nMap.GetEmax())
			sameEnergy = false;
	if (sameFOV)
		if (tMap.GetFovMin()!=nMap.GetFovMin() || tMap.GetFovMax()!=nMap.GetFovMax())
			sameFOV = false;
	if (sameCenter)
		if (tMap.GetMapCenterL()!=nMap.GetMapCenterL() || tMap.GetMapCenterB()!=nMap.GetMapCenterB())
			sameCenter = false;
	if (sameAlbedo)
		if (tMap.GetAlbedo()!=nMap.GetAlbedo())
			sameAlbedo = false;
	if (sameBinSize)
		if (tMap.GetYbin()!=nMap.GetYbin())
			sameBinSize = false;
	if (sameStep)
		if (m_expMaps[map].GetStep()!=m_expMaps[map+1].GetStep())
			sameStep = false;
	if (samePhCode)
		if (tMap.GetPhaseCode()!=nMap.GetPhaseCode())
			samePhCode = false;
	}


htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
htmlout << "<tr><td>Map</td><td>Counts</td><td>Date start</td><td>Date end</td><td>Energy</td><td>FOV</td>";
htmlout << "<td>Center</td><td>Albedo</td><td>BinSize</td><td>Step</td><td>PhCode</td></tr>" << endl;
for (int map=0; map<m_mapCount; ++map) {
	const AgileMap& m = m_ctsMaps[map];
	htmlout << "<tr><td>" << map+1 << "</td>";
	htmlout << "<td>" << m_counts[map] << "</td>";
	if (OpenConditionalBin(htmlout, sameDate, map, m_mapCount)) {
		WriteIsoDate(htmlout, m.GetStartDate());
		htmlout << "</td>";
		}
	if (OpenConditionalBin(htmlout, sameDate, map, m_mapCount)) {
		WriteIsoDate(htmlout, m.GetEndDate());
		htmlout << "</td>";
		}

	if (OpenConditionalBin(htmlout, sameEnergy, map, m_mapCount))
		htmlout << m.GetEmin() << ".." << m.GetEmax() << "</td>";
	if (OpenConditionalBin(htmlout, sameFOV, map, m_mapCount))
		htmlout << m.GetFovMin() << ".." << m.GetFovMax() << "</td>";

	if (OpenConditionalBin(htmlout, sameCenter, map, m_mapCount))
		htmlout <<m.GetMapCenterL() << ",&nbsp;" << m.GetMapCenterB() << "</td>";
	WriteConditionalBin(htmlout, m.GetAlbedo(), sameAlbedo, map, m_mapCount);
	WriteConditionalBin(htmlout, m.GetYbin(), sameBinSize, map, m_mapCount);
	WriteConditionalBin(htmlout, m_expMaps[map].GetStep(), sameStep, map, m_mapCount);
	WriteConditionalBin(htmlout, m.GetPhaseCode(), samePhCode, map, m_mapCount);
	htmlout << "</tr>" << endl;
/**
	htmlout << "</td><td>" << m.GetEmin()
		<< "</td><td>" << m.GetEmax()
		<< "</td><td>" << m.GetFovMin() << ".." << m.GetFovMax()
		<< "</td><td>" << m.GetMapCenterL() << ",&nbsp;" << m.GetMapCenterB()
		<< "</td><td>" << m.GetAlbedo()
		<< "</td><td>" << m.GetYbin()	/// binsize
		<< "</td><td>" << m_expMaps[map].GetStep()
		<< "</td><td>" << m.GetPhaseCode() << "</td></tr>" << endl;
*/
	}
htmlout << "</table>" << endl;
htmlout << "<br>" << endl;

/**
if (ExtCount()) {
	htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
	htmlout << "<tr><td>Map</td><td>Extended Source</td><td>Name</td><td>Flux</td><td>Index</td><td>FixFlag</td></tr>" << endl;
	for (int map=0; map<m_mapCount; ++map)
		for (int i=0; i<ExtCount(); ++i) {
			htmlout << "<tr><td>" << map+1 << "</td>";
			const AlikeDiffMap& extSrc = m_extSrcArr[ExtIndex(map, i)];
			char mapName[256];
			strcpy(mapName, extSrc.GetFileName());
			StripExtension(mapName, ".neb.gz", ".neb");
			htmlout << "<td>" << mapName << "</td>";
			/// MakeExtendedName(mapName, map, i);
			htmlout << "<td>" << mapName << "</td>";
			/// htmlout << "<td>" << m_initExtCoeffs[ExtIndex(map, i)] << "</td>";

			htmlout << "<td>" << m_extFlux[ExtIndex(map, i)] << "</td>";
			htmlout << "<td>" << m_extIndex[ExtIndex(map, i)] << "</td>";

			htmlout << "<td>" << (extSrc.GetFixflag()?"Free":"Fixed") << "</td></tr>" << endl;
			}
	htmlout << "</table>" << endl;
	htmlout << "<br>" << endl;
	}
*/

if (ExtCount()) {
	htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
	htmlout << "<tr><td>Extended Source</td><td>Barycentre</td><td>Flux</td><td>Index</td>";
	if (m_mapCount>1)
		htmlout << "<td>Files</td></tr>" << endl;
	else
		htmlout << "<td>File</td></tr>" << endl;
	for (int i=0; i<ExtCount(); ++i) {
		htmlout << "<tr><td>" << m_extData.Name(i) << "</td>";
		htmlout << "<td>" << m_extSrcArr[ExtIndex(0, i)].GetBaryL() << ", " << m_extSrcArr[ExtIndex(0, i)].GetBaryB() << "</td>";
		htmlout << "<td>" << fabs(m_extData.Flux(i));
		htmlout << (m_extSrcArr[ExtIndex(0, i)].GetFixflag()?" Free":" Fixed") << "</td>" << endl;
		htmlout << "<td>" << m_extData.Index(i) << "</td><td>";
		char mapName[256];
		const char* sep[2] = { "", ", " };
		for (int map=0; map<m_mapCount; ++map) {
			strcpy(mapName, m_extSrcArr[ExtIndex(map, i)].GetFileName());
			StripExtension(mapName, ".neb.gz", ".neb");
			htmlout << sep[bool(map)] << mapName;
			}
		htmlout << "</td></tr>" << endl;
		}
	htmlout << "</table>" << endl;
	htmlout << "<br>" << endl;
	}


if (SrcCount()) {
	htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
	htmlout << "<tr><td>Source</td><td>Flux</td><td>Index</td><td>L</td><td>B</td><td>sqrt(minTS)</td><td>FixFlag</td><td>ULCL</td><td>LOCL</td></tr>" << endl;
	for (int i=0; i<SrcCount(); ++i) {
		const SourceData& srcData = m_inSrcDataArr[i];
		htmlout << "<tr><td>" << srcData.label
					<< "</td><td>" << srcData.flux
					<< "</td><td>" << srcData.index
					<< "</td><td>" << srcData.srcL;
	
		if (m_srcLimits[i].lLim)
			htmlout << " in [" << m_srcLimits[i].lMin << ", " << m_srcLimits[i].lMax << "]";
	
		htmlout << "</td><td>" << srcData.srcB;
	
		if (m_srcLimits[i].bLim)
			htmlout << " in [" << m_srcLimits[i].bMin << ", " << m_srcLimits[i].bMax << "]";
	
		htmlout << "</td><td>" << sqrt(srcData.minTS) << "</td><td>" << srcData.fixflag << "</td><td>" << sqrt(m_sources[i].GetULCL()) << "</td><td>" << m_sources[i].GetLocCL() << "</td></tr>" << endl;
	
		}
	htmlout << "</table>" << endl;
	}


/// ///////////////////////////////////////////////////////

htmlout << "<h2>Output</h2>" << endl;

/// ///////////////////////////////////////////////////////


htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
htmlout << "<tr><td>DiffName</td><td>Coeff</td><td>Err</td><td>+Err</td><td>-Err</td>";

if (singleGal && singleIso) {
	AlikeDiffMap& gasMap = m_galSrc[0];
	AlikeDiffMap& isoMap = m_isoSrc[0];
	htmlout << "</tr>" << endl;
	htmlout << "<tr><td>Galactic</td><td>" << gasMap.GetCoeff() << "</td><td>" << gasMap.GetCoefferr()
			<< "</td><td>" << gasMap.GetCoeffhi() << "</td><td>" << gasMap.GetCoefflo() << "</td></tr>" << endl;
	htmlout << "<tr><td>Isotropic</td><td>" << isoMap.GetCoeff() << "</td><td>" << isoMap.GetCoefferr()
			<< "</td><td>" << isoMap.GetCoeffhi() << "</td><td>" << isoMap.GetCoefflo() << "</td></tr>" << endl;
	}
else {
	htmlout << "<td>DiffName</td><td>Coeff</td><td>Err</td><td>+Err</td><td>-Err</td></tr>" << endl;
	
	for (int map=0; map<m_mapCount; ++map) {
		AlikeDiffMap& gasMap = m_galSrc[map];
		AlikeDiffMap& isoMap = m_isoSrc[map];
		if (map==0) {
			const char* colSpan = "<td colspan=\"1\" rowspan=\"%d\" style=\"vertical-align: middle;\">";
			char galTD[128] = "<td>";
			char isoTD[128] = "<td>";
			if (singleGal)
				sprintf(galTD, colSpan, m_mapCount);
			else if (singleIso)
				sprintf(isoTD, colSpan, m_mapCount);

			htmlout << "<tr>";
			if (singleGal)
				htmlout << galTD << "Galactic" << "</td>";
			else
				htmlout << galTD << "Galactic_" << map+1 << "</td>";
			htmlout << galTD << gasMap.GetCoeff() << "</td>" << galTD << gasMap.GetCoefferr()
				<< "</td>" << galTD << gasMap.GetCoeffhi() << "</td>" << galTD << gasMap.GetCoefflo() << "</td>" << endl;

			if (singleIso)
				htmlout << isoTD << "Isotropic" << "</td>";
			else
				htmlout << isoTD << "Isotropic_" << map+1 << "</td>";
			htmlout << isoTD << isoMap.GetCoeff() << "</td>" << isoTD << isoMap.GetCoefferr()
				<< "</td>" << isoTD << isoMap.GetCoeffhi() << "</td>" << isoTD << isoMap.GetCoefflo() << "</td>" << endl;
			htmlout << "</tr>" << endl;
			}
		else {
			htmlout << "<tr>";
			if (!singleGal)
				htmlout << "<td>Galactic_" << map+1 << "</td><td>" << gasMap.GetCoeff() << "</td><td>" << gasMap.GetCoefferr()
				<< "</td><td>" << gasMap.GetCoeffhi() << "</td><td>" << gasMap.GetCoefflo() << "</td>" << endl;
			if (!singleIso)
				htmlout << "<td>Isotropic_" << map+1 << "</td><td>" << isoMap.GetCoeff() << "</td><td>" << isoMap.GetCoefferr()
				<< "</td><td>" << isoMap.GetCoeffhi() << "</td><td>" << isoMap.GetCoefflo() << "</td>" << endl;
			htmlout << "</tr>" << endl;
			}
		}
	}
htmlout << "</table>" << endl << "<br><br>" << endl;





if (ExtCount()) {
	htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
	htmlout << "<tr><td>Extended Source</td><td>Sqrt(TS)</td><td>Coeff</td><td>Err</td><td>+Err</td><td>-Err</td></tr>" << endl;
	for (int i=0; i<ExtCount(); ++i) {
		htmlout << "<tr><td>" << m_extData.Name(i) << "</td>";
		const AlikeExtMap& extSrc = m_extSrcArr[ExtIndex(0, i)];
		double BaryExposure = EvalExposure(extSrc.GetBaryL(), extSrc.GetBaryB(), m_expMaps[0]);
		htmlout << "<td>" << sqrt(extSrc.GetTS()) << "</td>";
		htmlout << "<td>" << extSrc.GetCoeff()/BaryExposure << "</td>";
		htmlout << "<td>" << extSrc.GetCoefferr()/BaryExposure << "</td>";
		htmlout << "<td>" << extSrc.GetCoeffhi()/BaryExposure << "</td>";
		htmlout << "<td>" << extSrc.GetCoefflo()/BaryExposure << "</td></tr>" << endl;
		}
	htmlout << "</table>" << endl;
	htmlout << "<br>" << endl;
	}

	
/// ///////////////////////////////////////////////////////


if (SrcCount()) {
	
	htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
	htmlout << "<tr><td>SrcName</td><td>sqrt(TS)</td><td>L_peak</td><td>B_peak</td><td>dist</td><td>Exp</td><td>Counts</td><td>Err</td><td>Flux</td><td>Err</td><td>Flux UL</td><td>Index</td><td>Err</td><td>L</td><td>B</td><td>distE</td><td>Radius</td><td>a</td><td>b</td><td>phi</td><td>gal</td><td>gal0</td><td>iso</td><td>iso0</td><td>fitinfo (cts fcn0 fcn1 edm0 edm1 iter0 iter1)</td></tr>" << endl;
	
	for (int i=0; i<m_srcCount; ++i) {
		const Ellipse& ellipse = m_sources[i].GetEllipse();
	/*
		double totalExposure = 0;
		for (int map=0; map<m_mapCount; ++map)
			totalExposure += m_sources[map*SrcCount()+i].GetExp();
	
		cerr << "**********************" << m_sources[i].GetExp() << endl;
		cerr << "**********************" << m_sources[i].GetCounts() << endl;
		cerr << "**********************" << m_sources[i].GetCountsErr() << endl;
	*/

		double exposure = GetTotalExposure(i);
		double dist = SphDistDeg(m_sources[i].GetSrcL(), m_sources[i].GetSrcB(), m_inSrcDataArr[i].srcL, m_inSrcDataArr[i].srcB);
		htmlout << "<tr><td>" << m_sources[i].GetLabel()
					<< "</td><td>" << sqrt(m_sources[i].GetTS())
					<< "</td><td>" << m_sources[i].GetSrcL()
					<< "</td><td>" << m_sources[i].GetSrcB()
		
					
					<< "</td><td>" << dist
/**
					<< "</td><td>" << m_sources[i].GetExp()
					<< "</td><td>" << m_sources[i].GetCounts()
					<< "</td><td>" << m_sources[i].GetCountsErr()
*/
					<< "</td><td>" << exposure
					<< "</td><td>" << exposure*m_sources[i].GetFlux()
					<< "</td><td>" << exposure*m_sources[i].GetFluxerr()

					<< "</td><td>" << m_sources[i].GetFlux()
					<< "</td><td>" << m_sources[i].GetFluxerr()
					<< "</td><td>" << m_sources[i].GetFluxul()
					<< "</td><td>" << m_sources[i].GetIndex()
		            << "</td><td>" << m_sources[i].GetIndexerr();
		
		
		if (ellipse.horAxis!=0 && ellipse.verAxis!=0) {
			double dist2 = SphDistDeg(ellipse.center.x, ellipse.center.y, m_inSrcDataArr[i].srcL, m_inSrcDataArr[i].srcB);
			htmlout << "</td><td>" << ellipse.center.x
					<< "</td><td>" << ellipse.center.y
					<< "</td><td>" << dist2
					<< "</td><td>" << m_sources[i].GetRadius()
					<< "</td><td>" << ellipse.horAxis
					<< "</td><td>" << ellipse.verAxis
			<< "</td><td>" << ellipse.attitude*RAD2DEG;
		} else {
			htmlout << "</td><td>-1</td><td>-1</td><td>-1</td><td>-1</td><td>-1</td><td>-1</td><td>-1";
		}
		
		if (m_inSrcDataArr[i].fixflag) {
			const char* sep[2] = {"", ","};
			htmlout << "</td><td>";
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPM(false, m, i, false);
			htmlout << " +/- ";
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPMErr(false, m, i, false);
			htmlout << "</td><td>";
			
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPM(false, m, i, true);
				htmlout << " +/- ";
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPMErr(false, m, i, true);
			htmlout << "</td><td>";
			
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPM(true, m, i, false);
				htmlout << " +/- ";
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPMErr(true, m, i, false);
			
			htmlout << "</td><td>";
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPM(true, m, i, true);
				htmlout << " +/- ";
			for (int m=0; m<m_mapCount; ++m)
				htmlout << sep[bool(m)] << GetFinalDPMErr(true, m, i, true);

		} else {
			htmlout << "</td><td>";
			
			htmlout << "</td><td>";
			
			htmlout << "</td><td>";
			
			htmlout << "</td><td>";
			
		}
		
		htmlout << "</td><td>";
		htmlout << m_fitInfo[i].counts;
		
		htmlout << " " << m_fitInfo[i].fcn0;
		htmlout << " " << m_fitInfo[i].fcn1;
		htmlout << " " << m_fitInfo[i].edm0;
		htmlout << " " << m_fitInfo[i].edm1;
		htmlout << " " << m_fitInfo[i].iter0;
		htmlout << " " << m_fitInfo[i].iter1;
		
			htmlout << "</td></tr>" << endl;
		
	}
	htmlout << "</table>" << endl;
	}


if (SrcCount() && m_mapCount>1) {
	bool multiChannel = false;
	double* eMinArr = new double[m_mapCount];
	double* eMaxArr = new double[m_mapCount];
	for (int map=0; map<m_mapCount; ++map) {
		eMinArr[map] = m_ctsMaps[map].GetEmin();
		eMaxArr[map] = m_ctsMaps[map].GetEmax();
		if (eMinArr[map]!=m_energyInf || eMaxArr[map]!=m_energySup)
			multiChannel = true;
		}
	if (multiChannel) {
		/// Sorting the energy channels
		for (int i=0; i<m_mapCount-1; ++i)
			for (int j=i+1; j<m_mapCount; ++j)
				if (eMinArr[i]>eMinArr[j] || (eMinArr[i]==eMinArr[j] && eMaxArr[i]>eMaxArr[j])) {
					double temp = eMinArr[i];
					eMinArr[i] = eMinArr[j];
					eMinArr[j] = temp;
					temp = eMaxArr[i];
					eMaxArr[i] = eMaxArr[j];
					eMaxArr[j] = temp;
					}
		/// Stripping redundant channels
		int chCount = m_mapCount;
		for (int i=0; i<chCount-1;)
			if (eMinArr[i]==eMinArr[i+1] && eMaxArr[i]==eMaxArr[i+1]) {
				for (int j=i+1; j<chCount-1; ++j) {
					eMinArr[j] = eMinArr[j+1];
					eMaxArr[j] = eMaxArr[j+1];
					}
				--chCount;
				}
			else
				++i;
		htmlout << endl << "<br><br>" << endl;
		htmlout << "<table border=1 cellpadding=2 cellspacing=0>" << endl;
		htmlout << "<tr><td>Source flux per channel</td>";
		for (int map=0; map<chCount; ++map)
			htmlout << "<td>" << eMinArr[map] << ".." << eMaxArr[map] << "</td>";
		htmlout << "</tr><tr>" << endl;
		for (int i=0; i<m_srcCount; ++i) {
			htmlout << "<tr><td>" << m_sources[i].GetLabel() << "</td>";
			double flux = m_sources[i].GetFlux();
			/// double index = m_sources[i].GetIndex();
			for (int map=0; map<chCount; ++map)
				htmlout << "<td>" << flux*m_sources[map*m_srcCount+i].GetNormFactor() << "</td>";
			htmlout << "</tr>" << endl;
			}
		htmlout << "</table>" << endl;
		}
	delete[] eMinArr;
	delete[] eMaxArr;
	}


htmlout << "</body></html>" << endl;
htmlout.close();
}
