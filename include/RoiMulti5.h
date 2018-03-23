//
// C++ Interface: %{MODULE}
//
// Description:
//
//
// Author: %{AUTHOR} <%{EMAIL}>, (C) %{YEAR}
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef _ROI_MULTI_
#define _ROI_MULTI_

#define PARNUM 6

/// #include <iomanip>
#include <fstream>

#include "TNamed.h"
#include "TVirtualFitter.h"
#include "TH1.h"
#include "TF1.h"
#include "TMinuit.h"

#include "AlikeData5.h"
#include "AlikePsf.h"
#include "Ellipse.h"

/* 03/2018  -> ExpRatioEvaluator dependency */
#include "ExpRatioEvaluator.h"


#define ROIMULTI_VERSION "RoiMulti 2.0"
#define ROIMULTI_DATE "2012-06-21"




Double_t EvalExposure(Double_t srcl, Double_t srcb, const AgileMap& expmap);



class AlikeDiffMap: public AgileMap, public AlikeNorm
{
public:	/// Contruction and Copy
	AlikeDiffMap(): AgileMap(), AlikeNorm(), m_index(2.0), m_coeff(1), m_coefferr(0), m_coefflo(0), m_coeffhi(0), m_fixflag(1) { UpdateNorm(); }
	AlikeDiffMap(const AlikeDiffMap &another): AgileMap(), AlikeNorm() { Copy(another, false); }
	~AlikeDiffMap() {}

	void BuildFrom(const AgileMap &map, double eInf, double eSup, double index=2.0);

	AlikeDiffMap& operator=(const AlikeDiffMap& another) { if (this!=&another) Copy(another, true); return *this; }

public:	/// Get Functions
	Double_t GetIndex() const { return m_index; }
	Double_t GetCoeff() const { return m_coeff; }
	Double_t GetCoefferr() const { return m_coefferr; }
	Double_t GetCoefflo() const { return m_coefflo; }
	Double_t GetCoeffhi() const { return m_coeffhi; }
	Int_t    GetFixflag() const { return m_fixflag; }


public:	/// Set Functions
	void SetIndex(Double_t index) { m_index = index>0?index:-index; UpdateNorm(); }
	void SetCoeff(Double_t coeff) { m_coeff = coeff; }
	void SetCoefferr(Double_t coefferr) { m_coefferr = coefferr; }
	void SetCoefflo(Double_t coefflo) { m_coefflo = coefflo; }
	void SetCoeffhi(Double_t coeffhi) { m_coeffhi = coeffhi; }
	void SetFixflag(Int_t fixflag) { m_fixflag = fixflag; }

	void Expose(const AgileMap& expMap); /// multiply by the exposition

private:	/// Data
	Double_t m_index;
	Double_t m_coeff;
	Double_t m_coefferr;
	Double_t m_coefflo;
	Double_t m_coeffhi;
	FixFlag  m_fixflag;
	void UpdateNorm() { AlikeNorm::UpdateNorm(GetEmin(), GetEmax(), m_index); }


protected:	/// Making a copy
	void Copy(const AlikeDiffMap &another, bool all);

private:	/// Internal Operations
	void UpdNorm();
};


class AlikeExtMap: public AlikeDiffMap
{
public:
	AlikeExtMap(): AlikeDiffMap(), m_TS(0), m_baryl(0), m_baryb(0) {}
	AlikeExtMap(const AlikeExtMap &another): AlikeDiffMap(another) { Copy(another, false); }

	void BuildFrom(const AgileMap& map, double eInf, double eSup, double index=2.0)
		{ AlikeDiffMap::BuildFrom(map, eInf, eSup, index); m_TS=0; EvalBarycenter(); }
	AlikeExtMap& operator=(const AlikeExtMap& another) { if (this!=&another) Copy(another, true); return *this; }

	bool CompatibleWith(const AgileMap& map) const { return SameParamsAs(map); }
	double DegDistance(double l, double b) const { return SphDistDeg(GetBaryL(), GetBaryB(), l, b); }

	double GetTS() const { return m_TS; }
	void SetTS(double ts) ;
	void SetLikelihood(Double_t lik) { m_likelihood = lik; };
	double GetLikelihood() { return m_likelihood; };
	void PrintSqrtTS() const;
	double GetBaryL() const { return m_baryl; }
	double GetBaryB() const { return m_baryb; }

private:	/// Data
	double   m_TS;
	double	 m_likelihood;
	double   m_baryl;
	double   m_baryb;
	void EvalBarycenter();

protected:	/// Making a copy
	void Copy(const AlikeExtMap &another, bool all);
};


class AlikeSourceMap: public AlikePsfSource
{
public: /// Construction
	AlikeSourceMap(): AlikePsfSource(),
		m_flux(0), m_fluxerr(0), m_fluxlo(0), m_fluxhi(0), m_fluxul(0),
		m_idxErr(0), m_par2Err(0), m_par3Err(0), m_exposure(0), m_ts(0), m_minTS(3),
		m_ulcl(0.95), m_loccl(0.95), m_fixflag(0), m_forceposfree(0), /// m_covar(),
		m_polygon(), m_ellipse(), m_radius(0), m_label()
        {
            for(int i=0; i<7; i++) {
                m_status[i] = -100;
                m_cts[i] = -100;
            }
        }

	~AlikeSourceMap() {}

	void Set(
		const AlikePsfTables* psfTab,
		const AgileMap &map,
		double eInf, double eSup, double theta,
		double srcL, double srcB, double index,
		int typefun, double par2, double par3,
		Double_t flux,
		Double_t exposure)
		{
		AlikePsfSource::Set(psfTab, map, eInf, eSup, theta, srcL, srcB, index, typefun, par2, par3);
		m_flux = flux;
		m_exposure = exposure;
		}

public:	/// Data Access
	Double_t GetFlux() const { return m_flux; }
	Double_t GetFluxerr() const { return m_fluxerr; }
	Double_t GetFluxlo() const { return m_fluxlo; }
	Double_t GetFluxhi() const { return m_fluxhi; }
	Double_t GetFluxul() const { return m_fluxul; }
	void SetFlux(Double_t flux) { m_flux = flux; }
	void SetFluxerr(Double_t fluxerr) { m_fluxerr = fluxerr; }
	void SetFluxlo(Double_t fluxlo) { m_fluxlo = fluxlo; }
	void SetFluxhi(Double_t fluxhi) { m_fluxhi = fluxhi; }
	void SetFluxul(Double_t fluxul) { m_fluxul = fluxul; }

/**
	Double_t GetCounts() const { return m_flux*m_exposure; }
	Double_t GetCountsErr() const { return m_fluxerr*m_exposure; }
	Double_t GetCountslo() const { return m_fluxlo*m_exposure; }
	Double_t GetCountshi() const { return m_fluxhi*m_exposure; }
	Double_t GetCountsul() const { return m_fluxul*m_exposure; }
*/
	void SetIndexerr(double idxErr) { m_idxErr = idxErr; }
	double GetIndexerr() const { return m_idxErr; }

	void SetPar2err(double idxErr) { m_par2Err = idxErr; }
	double GetPar2err() const { return m_par2Err; }

	void SetPar3err(double idxErr) { m_par3Err = idxErr; }
	double GetPar3err() const { return m_par3Err; }

	Double_t GetExp() const { return m_exposure; }
	void SetExp(Double_t exposure) { m_exposure = exposure; }

	Double_t GetMinTS() const { return m_minTS; }
	void SetMinTS(Double_t minTS) { m_minTS = minTS; }

	Int_t GetFixflag() const { return m_fixflag; }
	void SetFixflag(Int_t fixflag) { m_fixflag = fixflag; if(fixflag & ForcePosFree) m_forceposfree = 1; }
	Int_t   GetForcePosFree() const { return m_forceposfree; }

	Double_t GetTS() const { return m_ts; }
	void SetTS(Double_t inTS) ;

	void SetLikelihood(Double_t lik) { m_likelihood = lik; };
	Double_t GetLikelihood() { return m_likelihood; };

	Double_t GetULCL() const { return m_ulcl; }
	Double_t GetLocCL() const { return m_loccl; }

	void SetULCL(Double_t ulcl) ;
	void SetLocCL(Double_t loccl) ;
/**
	const TMatrixD& GetCovar() const { return m_covar; }
	void SetCovar(const TMatrixD& covar) { m_covar.ResizeTo(covar); m_covar = covar; }
*/
	const Polygon& GetPoly() const { return m_polygon; }
	const Ellipse& GetEllipse() const { return m_ellipse; }
	double GetRadius() const { return m_radius; }
	void SetContour(const Polygon& polygon, const Ellipse& ellipse, double radius)
				{ m_polygon=polygon; m_ellipse=ellipse; m_radius=radius; }
	void WritePolygon(string filesuffix) const;
	void WriteEllipse(string filesuffix) const;

	void PrintSqrtTS() const;
	bool PrintAbphi() const;	/// return true if something was printed

	string GetLabel() const { return m_label; }
	string SetLabel(string label) { m_label = label; return m_label; }

    void SetStatus(int step, int status) { m_status[step] = status; }
    int GetStatus(int step) { return m_status[step]; }

    void SetCts(int step, int cts) { m_cts[step] = cts; }
    int GetCts(int step) { return m_cts[step]; }

	double GetSpectraCorrectionFactor(bool fluxcorrection);

private: /// Data
	Double_t m_flux;
	Double_t m_fluxerr;
	Double_t m_fluxlo;
	Double_t m_fluxhi;
	Double_t m_fluxul;

	double m_idxErr;
	double m_par2Err;
	double m_par3Err;

	Double_t m_exposure;

	Double_t m_ts;
	Double_t m_minTS;
	Double_t m_likelihood;

	Double_t m_ulcl;
	Double_t m_loccl;
	FixFlag  m_fixflag;
	Int_t m_forceposfree;
	/// TMatrixD m_covar;

	Polygon  m_polygon;	/// Source coutour if available, empty otherwise
	Ellipse  m_ellipse;	/// Ellipse centered in m_polygon barycenter and best fitting it
	double   m_radius;	/// Radius of a circle with the same area as m_polygon
	string   m_label;

	int m_status[7]; // Fit status for each one of the 7 analysis steps
	int m_cts[7];    // Last counts of cts in move for the source, for each one of the 7 analysis steps

private: /// No implementation
	AlikeSourceMap(const AlikeSourceMap&);
	AlikeSourceMap& operator=(const AlikeSourceMap&);

};



/// Structure storing the limits in L and B parameters for the sources
struct SourceLimits
{
double lMin;
double lMax;
double bMin;
double bMax;
bool   lLim;	/// L is limited in [lMin..lMax]
bool   bLim;	/// B is limited in [bMin..bMax]
};


/**
void InitSimulationEngine();
AgileMap* NewSimulationArray(const AgileMap& model, int count);
*/

typedef int DiffMode;
const DiffMode DiffFixed      = 0;
const DiffMode DiffDefault    = 1;
const DiffMode DiffFree       = 2;
const DiffMode DiffTogether   = 3;




struct FitInfo
{
	//amin or fcn : chisquare
	//edm : estimated distance to minimum
	//errdef
	//nvpar : number of variable parameters
	//nparx : total number of parameters
	FitInfo(): fcn0(0), fcn1(0), edm0(0), edm1(0), iter0(0), iter1(0), counts(0), fitstatus0(0), fitstatus1(0), nvpar0(0), nvpar1(0), nparx0(0), nparx1(0) {}
	double fcn0;
	double fcn1;
	double edm0;
	double edm1;
	int    iter0;
	int    iter1;
	int fitstatus0;
	int fitstatus1;
	int nvpar0;
	int nvpar1;
	int nparx0;
	int nparx1;
	int    counts;
};




/// Statistical settings

unsigned GetSeed();
unsigned SetSeed(unsigned seed);





/// ////////////////////////////////////////////////////
///
/// class FitParameters
///
/// Class for remapping the model parameters indices
///
/// This class stores an array of values, some constant and some variable
/// It gives access either to the entire array or to the variable subset




class FitParameters
{
public:
	FitParameters(): m_params(0), m_vars(0), m_count(0), m_map(0), m_mapCount(0) {}
	~FitParameters() { delete[] m_params; delete[] m_vars; delete[] m_map; }

/// Setting up the values (complete array)
	void Init(int count, bool allVariables = true);

	void SetValue(int i, double value) { m_params[i]=value; }
	void SetValue(int i, double value, bool makeVariable);

	void FixValue(int i, double value) { SetValue(i, value, false); }
	void FreeValue(int i, double value) { SetValue(i, value, true); }

	double& operator[](int i) { return m_params[i]; }

/// Setting up the values (variable subset)
	void UpdateValues(const double* varVals);	/// given a shorter array with the variable values only
	void SetVarValue(int i, double varVal) { m_params[m_map[i]] = varVal; } /// index of a variable value

/// Getting the values (complete array)
	int Count() const { return m_count; }
	bool IsVariable(int i) const { return m_vars[i]; }
	double operator[](int i) const { return m_params[i]; }

	int VariableCount() const { return m_mapCount; }
	int VariableIndex(int i) const { return m_map[i]; }
	double GetVarValue(int i) const { return m_params[m_map[i]]; }
	/// Given the incomplete array varVals, fill the full size array allVals
	void FillParArray(const double* varVals, double* allVals) const;

/// Setting and getting the TF1 parameters
	void SetModelValues(TF1& model) const;
	void GetModelValues(const TF1& model);

private:
	double* m_params; 	/// Values
	bool*   m_vars;		/// Variable flags
	int     m_count;		/// Number of values
	int*    m_map;			/// Map of the variable values
	int     m_mapCount;	/// Number ofvariable values (length of m_map array)
	void UpdateMap();
};




/// Same as BUILD21
#define FLUX_SCALE_FACTOR 1000000.0
#define FLUX_UPPER_BOUND 10000000.0


/** other factors tried with little success:

#define FLUX_SCALE_FACTOR 100000.0
#define FLUX_UPPER_BOUND 10.0

#define FLUX_SCALE_FACTOR 1.0
#define FLUX_UPPER_BOUND 10.0

#define FLUX_SCALE_FACTOR 1.0
#define FLUX_UPPER_BOUND 100.0

#define FLUX_SCALE_FACTOR 1.0
#define FLUX_UPPER_BOUND 100000.0
*/


class RoiMulti
{
static RoiMulti* sm_singleton;	/// Only one instance may exist
static int       s_fitterIterations;
static unsigned* s_fitterChanges;
public:	/// Construction
	RoiMulti();
	~RoiMulti();

public:	/// Getting the singleton object
	static RoiMulti& Get() { return *sm_singleton; }

public:	/// Main operations
	bool SetPsf(const char* psfFileName, const char* raeffFileName, const char* edpFileName);
	bool SetMaps(const MapData& mapData, int galMode=DiffDefault, int isoMode=DiffDefault);
	bool SetMinimizer(const char* minimizertype, const char* minimizeralg, int minimizerdefstrategy, double deftol);
	bool SetCorrections(int galmode2, int galmode2fit, int isomode2, int isomode2fit, int edpcorrection, int fluxcorrection) { m_galmode2 = galmode2; m_galmode2fit = galmode2fit; m_isomode2 = isomode2; m_isomode2fit = isomode2fit; m_edpcorrection = edpcorrection; m_fluxcorrection = fluxcorrection;};
	/// Add the Extended Sources for analysis
	bool SetExtendedSources(const ExtData& extData);

	/// Analysis
	bool SetLogfile(const char* fileName);

	void CloseLogFile();
	int DoFit(
		const SourceDataArray& srcDataArr,
		double ranal,
		double ulcl,
		double loccl,
		int chatter = 1,
		const char* sourceCheckTS = 0,
		double minSourceTS = 0);
	static int GetFitterIterations() { return s_fitterIterations; }
	void MakeTSArrays(double*& tsVarArr, double*& tsConstArr);


	/// Simulation
	/// AgileMap MakeModel(const SourceDataArray& srcDataArr);
	AgileMap* NewModelArray(const SourceDataArray& srcDataArr);
	AgileMap* NewSimulationArray(const SourceDataArray& srcDataArr);
	bool WriteSimulatedCounts(const SourceDataArray& srcDataArr, const char* outfilename);
	bool WriteSimulatedCounts(const SourceDataArray& srcDataArr, const MapList& fileNames);
	bool WriteModelCounts(const SourceDataArray& srcDataArr, const char* outfilename);
	bool WriteModelCounts(const SourceDataArray& srcDataArr, const MapList& fileNames);

public:	/// Collecting the results
	void Write(const char* fileName, bool skipFitInfo=true) const;
	void WriteSources(const char* fileName, bool expratioevaluation, bool isExpMapNormalized=false, double minThreshold=0, double maxThreshold=0, int squareSize=0, bool skipFixed=false, bool skipEllipses=false) const;
	void LogSources(const char* fileName, int iterNum, AgileMap* simArr, int simArraySize) const;
	void WriteHtml(const char* fileName, bool expratioevaluation, bool isExpMapNormalized=false, double minThreshold=0, double maxThreshold=0, int squareSize=0, const char* suffix=".html") const;
	SourceDataArray GetFitData() const { return m_outSrcDataArr; }

	/// Diffuse components
	int MapCount() const { return m_mapCount; }
	/// int GalCount() const { return m_galMode==3?1:m_mapCount; }
	/// int IsoCount() const { return m_isoMode==3?1:m_mapCount; }
	int ExtCount() const { return m_extCount; }
	const AlikeDiffMap& GetGalactic(int index) const { return m_galSrc[index%m_mapCount]; }
	const AlikeDiffMap& GetIsotropic(int index) const { return m_isoSrc[index%m_mapCount]; }
	const AlikeExtMap& GetExtended(int map, int index) const { return m_extSrcArr[ExtIndex(map, index)]; }

	/// Sources
	int SrcCount() const { return m_srcCount; }
	const AlikeSourceMap& GetSource(int source) const { return m_sources[source%m_srcCount]; }
	double GetTotalExposure(int source) const;
	double GetTotalExposureSpectraCorrected(int source) const;
	Double_t GetSourceTS(const char* label) const;

private:	/// Internal operations
	void SetSources(const SourceDataArray& srcDataArr, double ranal, double ulcl, double loccl);
	Int_t DoFit(const char* fitOpt1, const char* fitOpt2, const char* sourceCheckTS, double minSourceTS);
	void FixDistantSources(double L, double B, int source, bool extended);
	void FitExtendedSources(int source, const char* fitOpt, bool* diffParsStored);
	void FitFlux(int activeSource, const char* fitOpt);
	void FitPositionIndex(int activeSource, const char* fitOpt);
	bool LoopExt(const char* fitOpt);
	void Loop1(const char* fitOpt);
	void Loop2(const char* fitOpt);
	void MakeCovarMat(int source);
	void MakeEllipse(int source);

	/// Cleaning memory
	void CleanMaps();
	void CleanSources();

	/// Parameters access

	/// map in [0..m_mapCount-1]
	int GalCoeffPar(int map) const { return m_galMode==3?0:map; }
	int IsoCoeffPar(int map) const { return m_galParCount+(m_isoMode==3?0:map); }

	/// int ExtCoeffPar(int map, int index) const { return m_diffParCount+ExtIndex(map,index); }
	int ExtCoeffPar(int index) const { return m_diffParCount+index; }

	/// src in [0..m_srcCount-1]
	int SrcFluxPar(int src) const { return m_sourceParOffset + src*6; }
	int SrcIdxPar(int src) const  { return m_sourceParOffset + src*6 + 1; }
	int SrcLPar(int src) const    { return m_sourceParOffset + src*6 + 2; }
	int SrcBPar(int src) const    { return m_sourceParOffset + src*6 + 3; }
	int SrcPar2Par(int src) const  { return m_sourceParOffset + src*6 + 4; }
	int SrcPar3Par(int src) const  { return m_sourceParOffset + src*6 + 5; }

	/// Fitting
	Double_t FitFunction(Double_t *x, Double_t *par);
	static Double_t SingleFitFunction(Double_t *x, Double_t *par) { return sm_singleton->FitFunction(x, par); }
	Int_t Fit(const char* opt, int source, bool zero, int flags, Double_t* amin=0);
	double Likelihood(double* data=0);

	/// Setting up the Model
	void InitParams();
	void Move(double lCenter, double bCenter);

	/// Diffuse parameters
    void SetGalPar(int map, double val) { m_model.SetParameter(GalCoeffPar(map), val); }
	void SetIsoPar(int map, double val) { m_model.SetParameter(IsoCoeffPar(map), val); }
	void SetGalPar(int map) { SetGalPar(map, m_galSrc[map].GetCoeff()); }
	void SetIsoPar(int map) { SetIsoPar(map, m_isoSrc[map].GetCoeff()); }

	void FixGalPar(int map) { m_model.FixParameter(GalCoeffPar(map), m_galSrc[map].GetCoeff()); }
	void FixIsoPar(int map) { m_model.FixParameter(IsoCoeffPar(map), m_isoSrc[map].GetCoeff()); }

	void ReleaseGalPar(int map) { m_model.SetParLimits(GalCoeffPar(map), m_galLimitMin, m_galLimitMax); }
	void ReleaseIsoPar(int map) { m_model.SetParLimits(IsoCoeffPar(map), m_isoLimitMin, m_isoLimitMax); }

	/// Extended Sources
    void NullifyExt(int source) { m_model.FixParameter(ExtCoeffPar(source), 0); }
	void SetExtPar(int source, double val) { m_model.SetParameter(ExtCoeffPar(source), val); }
	void SetExtPar(int source) { SetExtPar(source, m_extSrcArr[ExtIndex(0, source)].GetCoeff()); }
	void FixExtPar(int source) { m_model.FixParameter(ExtCoeffPar(source), m_extSrcArr[ExtIndex(0, source)].GetCoeff()); }
	void ReleaseExtPar(int source) { m_model.SetParLimits(ExtCoeffPar(source), m_extLimitMin, m_extLimitMax); }
	void GetExtPar(int source);

	/// Source parameters
	void Nullify(int source);
	void SetSrcPars(int source);

	/// Flux parameter
    void ReleaseSrcFlux(int source) { m_model.SetParLimits(SrcFluxPar(source), m_fluxLimitMin, m_fluxLimitMax); }
	void SetSrcFlux(int source, double flux) { m_model.SetParameter(SrcFluxPar(source), flux*m_fluxScaleFactor); }
	void SetSrcFlux(int source) { SetSrcFlux(source, m_sources[source].GetFlux()); }
	void FixSrcFlux(int source, double flux) { m_model.FixParameter(SrcFluxPar(source), flux*m_fluxScaleFactor); }
	void FixSrcFlux(int source) { FixSrcFlux(source, m_sources[source].GetFlux()); }

	/// Position parameters
	void ReleaseSrcPos(int source);
	void SetSrcPos(int source) { SetSrcPos(source, m_sources[source].GetSrcL(), m_sources[source].GetSrcB()); }
	void SetSrcPos(int source, double l, double b);
	void FixSrcPos(int source) { FixSrcPos(source, m_sources[source].GetSrcL(), m_sources[source].GetSrcB()); }
	void FixSrcPos(int source, double l, double b);

	/// Spectral index parameter
    void ReleaseSrcIndex(int source) { m_model.SetParLimits(SrcIdxPar(source), m_indexLimitMin, m_indexLimitMax); }
	void SetSrcIndex(int source, double index) { m_model.SetParameter(SrcIdxPar(source), index); }
	void SetSrcIndex(int source) { SetSrcIndex(source, m_sources[source].GetIndex()); }
	void FixSrcIndex(int source, double index) { m_model.FixParameter(SrcIdxPar(source), index); }
	void FixSrcIndex(int source) { FixSrcIndex(source, m_sources[source].GetIndex()); }

	/// Par2 parameter
	void ReleaseSrcPar2(int source) { m_model.SetParLimits(SrcPar2Par(source), m_par2LimitMin, m_par2LimitMax); }
	void SetSrcPar2(int source, double par2) { m_model.SetParameter(SrcPar2Par(source), par2); }
	void SetSrcPar2(int source) { SetSrcPar2(source, m_sources[source].GetPar2()); }
	void FixSrcPar2(int source, double par2) { m_model.FixParameter(SrcPar2Par(source), par2); }
	void FixSrcPar2(int source) { FixSrcPar2(source, m_sources[source].GetPar2()); }

	/// Par3 parameter
	void ReleaseSrcPar3(int source) { m_model.SetParLimits(SrcPar3Par(source), m_par3LimitMin, m_par3LimitMax); }
	void SetSrcPar3(int source, double par3) { m_model.SetParameter(SrcPar3Par(source), par3); }
	void SetSrcPar3(int source) { SetSrcPar3(source, m_sources[source].GetPar3()); }
	void FixSrcPar3(int source, double par3) { m_model.FixParameter(SrcPar3Par(source), par3); }
	void FixSrcPar3(int source) { FixSrcPar3(source, m_sources[source].GetPar3()); }

	/// Reading the Model and setting up diff and source components
	double PeekSrcFluxPar(int source) const { return m_model.GetParameter(SrcFluxPar(source))*m_fluxScaleFactorMulInv; }
	void GetSrcPars(int source);
	void GetFluxErrors(int source);

	void FitIndex(int source);
	void GetSrcUL(int source);

	void GetDiffPars();
	void GetDiffErrs();
	void PrintDiffData();

	double EvalFitFunction(double* params, double* data=0);

	// exp ratio evaluation - 03/2018
	double ExpRatioEvaluation(AgileMap& exp, double l, double b, bool isExpMapNormalized, double minThreshold, double maxThreshold, int squareSize) const;


private:	/// Data
	TH1D   m_countsHist;
	TF1    m_model;
	/// FitParameters m_fitParams;

	double m_ranal;
	double m_ulcl;
	double m_loccl;

	double m_lCenter;
	double m_bCenter;

	/// Psf
	AlikePsfTables m_psfTab;

	/// Maps
	int       m_mapCount;
	AgileMap* m_ctsMaps;
	AgileMap* m_expMaps;
	AgileMap* m_gasMaps;
	bool      m_alignedMaps;
	int*      m_counts;
	double*   m_thetaArr;
	double*   m_galCoeffs;
	double*   m_isoCoeffs;
	double    m_energyInf;
	double    m_energySup;

	/// Corrections
	int m_galmode2;
	int m_galmode2fit;
	int m_isomode2;
	int m_isomode2fit;
	int m_edpcorrection;
	int m_fluxcorrection;

	//Fitter
	int m_minimizerdefstrategy;

	/// Sources
	SourceDataArray m_inSrcDataArr;	/// Copy of the original input data
	SourceDataArray m_outSrcDataArr;	/// Input data modified after fitting

	/// Diff components
	DiffMode        m_galMode;
	DiffMode        m_isoMode;
	AlikeDiffMap*   m_galSrc;
	AlikeDiffMap*   m_isoSrc;
	int             m_galParCount;
	int             m_isoParCount;
	int             m_diffParCount; /// m_galParCount+m_isoParCount
	bool SingleGalPar() const;
	bool SingleIsoPar() const;

	void PrintDiffParams() const;

	Mat3D m_diffPars[2];
	Mat3D m_diffErrs[2];

	double GetDPM(bool iso, int par, int source, bool zero) const { return m_diffPars[iso](par, source, int(zero)); }
	void SetDPM(bool iso, int par, int source, bool zero, double val) { m_diffPars[iso](par, source, int(zero))=val; }
	double GetDPMErr(bool iso, int par, int source, bool zero) const { return m_diffErrs[iso](par, source, int(zero)); }
	void SetDPMErr(bool iso, int par, int source, bool zero, double err) { m_diffErrs[iso](par, source, int(zero))=err; }

	double GetFinalDPM(bool iso, int par, int source, bool zero) const; /// Multiplied by m_galCoeffs[par] or m_isoCoeffs[par] when the case
	double GetFinalDPMErr(bool iso, int par, int source, bool zero) const; /// Multiplied by m_galCoeffs[par] or m_isoCoeffs[par] when the case

    void ResetFitStatus() { m_status = -1; }
    void ResetFitCts() { m_cts = -1; }
    // Stop status assignment on errors (with last cts counts) for each fit loop.
    void SetFitStatus(int status) { if(m_status <= 0) { m_status = status; m_cts = m_ctsMove; } }

	/// Extended sources
	ExtData         m_extData;
	int             m_extCount;
	int             ExtIndex(int map, int index) const { return map*m_extCount+index; }
	AlikeExtMap*    m_extSrcArr; /// Length of this is m_mapCount*m_extCount
	double*         m_extFlux;
	double*         m_extIndex;
	/// double*         m_initExtCoeffs;

	/// Source components
	int             m_srcCount;
	AlikeSourceMap* m_sources;
	SourceLimits*   m_srcLimits;
	FitInfo*        m_fitInfo;
	int             m_sourceParOffset; /// m_diffParCount+m_mapCount*m_extCount

	std::ofstream*       m_logFile;
/**
	void Log(const char* text);
	ostream& Logfile() { return m_logFile ? *m_logFile : std::clog; }
*/

	/// Map references
	int*            m_binList;
	int*            m_binAddr;
	int             m_binCount;

	/// Tuning the flux parameter
	double m_fluxScaleFactor;
	double m_fluxScaleFactorMulInv;
	double m_fluxUpperBound;

    int m_status;
    int m_cts;
    int m_ctsMove;

    double m_galLimitMin;
    double m_galLimitMax;
    double m_isoLimitMin;
    double m_isoLimitMax;
    double m_extLimitMin;
    double m_extLimitMax;
    double m_fluxLimitMin;
    double m_fluxLimitMax;
    double m_indexLimitMin;
    double m_indexLimitMax;
	double m_par2LimitMin;
	double m_par2LimitMax;
	double m_par3LimitMin;
	double m_par3LimitMax;

private:		/// Forbidden functions
	RoiMulti(const RoiMulti& another);					/// Singleton class, copy constructor not allowed
	RoiMulti& operator=(const RoiMulti& another);	/// Singleton class, copy operator not allowed
};

#endif
