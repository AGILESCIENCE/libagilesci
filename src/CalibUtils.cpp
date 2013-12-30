////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       Scientific pipeline I/O routine
//       CalibUtils.cpp
//       Release: V1.0 -  14/Feb/2012
//       Contributors: 
//       Author: Andrew Chen, Tomaso Contessi (IASF-Milano), Andrea Bulgarelli (IASF-Bologna)
//
// INPUT
//       TBD
//
// OUTPUT
//       TBD
//
//
// FILE HISTORY
//       28/May/2012
//             First release: V1.0
//             Author: Andrew Chen, Tomaso Contessi (IASF-Milano), Andrea Bulgarelli (IASF-Bologna)
//             Based on 3 classes previously defined in AlikeLib.h
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>
#include <iostream>

/// #include "fitsio.h"
#include "FitsUtils.h"

#include "CalibUtils.h"


using namespace std;


typedef ArrayOf<1,long> VecL;


static void ReadFloatArray(FitsFile& file, const char* colName, VecF& floatArr)
{
if (file.Status())
	return;
int colnum = 0;
fits_get_colnum(file, CASEINSEN, const_cast<char*>(colName), &colnum, file);
int typecode = 0;
long nrows = 0;
long width = 0;
fits_get_eqcoltype(file, colnum, &typecode, &nrows, &width, file);
floatArr.ResizeTo(nrows);
fits_read_col(file, TFLOAT, colnum, 1, 1, nrows, NULL, floatArr/*.Buffer()*/, NULL, file);
}


static bool GetImageSize(FitsFile& file, VecL& naxes)
{
int naxis = 0;
fits_get_img_dim(file, &naxis, file);
if (file.Ok()) {
	naxes.ResizeTo(naxis);
	fits_get_img_size(file, naxis, naxes/*.Buffer()*/, file);
	}
return file.Ok();
}


/**
static bool ReshapeImage(FitsFile& file, MatF& image)
{
MatrixOf<long> naxes;
if (GetImageSize(file, naxes)) {
	MatrixOf<int> shape(naxis+1);
	for (int i=0; i<naxis; ++i)
		shape[i] = naxes[naxis-i-1];
	shape[naxis] = 0;
	image.ReshapeTo(shape.Buffer());
	return true;
	}
return false;
}
*/


/**
static void ReadFloatArray(fitsfile* file, const char* colName, int* status, VecF& floatArr)
{
if (*status)
	return;
long nrows = 0;
int colnum = 0;
int typecode = 0;
long width = 0;

fits_get_colnum(file, CASEINSEN, const_cast<char*>(colName), &colnum, status);
fits_get_eqcoltype(file, colnum, &typecode, &nrows, &width, status);
floatArr.ResizeTo(nrows);
fits_read_col(file, TFLOAT, colnum, 1, 1, nrows, NULL, floatArr.Buffer(), NULL, status);
}
*/


/// /////////////////////////////////////////////////////////////////////
///
/// class AeffGrid
///
///
///
///

int AeffGrid::Read(const char* aeffFileName)
{
int status = 0;
FitsFile reffFile;
if (!reffFile.Open(aeffFileName)) {
	cerr << "Error " << status << " opening file" << aeffFileName << endl;
	return status;
	}

reffFile.MoveAbsHDU(1);
VecL naxes;
if (GetImageSize(reffFile, naxes)) {
	int anynul = 0;
	m_aeffgrid.ResizeTo(naxes[2], naxes[1], naxes[0]);
	fits_read_3d_flt(reffFile, 0, 0, naxes[2], naxes[1], naxes[2], naxes[1], naxes[0], m_aeffgrid/*.Buffer()*/, &anynul, reffFile);
	
	reffFile.MoveAbsHDU(3);
	ReadFloatArray(reffFile, "ENERGY", m_energy);
	ReadFloatArray(reffFile, "POLAR_ANGLE", m_theta);
	ReadFloatArray(reffFile, "AZIMUTH_ANGLE", m_phi);
	}
return reffFile.Status();
}

/**
int AeffGrid::Read(const char* aeffFileName)
{
int status = 0;
fitsfile* reffFile;
if (fits_open_file(&reffFile, aeffFileName, READONLY, &status)) {
	cerr << "Error " << status << " opening file" << aeffFileName << endl;
	return status;
	}

int  naxis = 0;
/// int  bitpix = 0; zzz unused
long naxes[3];

fits_movabs_hdu(reffFile, 1, NULL, &status);
/// fits_get_img_type(reffFile, &bitpix, &status); zzz unused
fits_get_img_dim(reffFile, &naxis, &status);
fits_get_img_size(reffFile, naxis, naxes, &status);

if (status==0) {
	m_aeffgrid.ReshapeTo(3, naxes[2], naxes[1], naxes[0]);
	
	int anynul = 0;
	fits_read_3d_flt(reffFile, 0, 0, naxes[2], naxes[1], naxes[2], naxes[1], naxes[0], m_aeffgrid.Buffer(), &anynul, &status);
	
	fits_movabs_hdu(reffFile, 3, NULL, &status);
	
	ReadFloatArray(reffFile, "ENERGY", &status, m_energy);
	ReadFloatArray(reffFile, "POLAR_ANGLE", &status, m_theta);
	ReadFloatArray(reffFile, "AZIMUTH_ANGLE", &status, m_phi);
	}
int closeStatus = 0;
fits_close_file(reffFile, &closeStatus);

return status;
}
*/


static void FindInterval(const VecF& arr, float value, int* minIndex, int* maxIndex)
{

*minIndex = arr.LeftIndex(value);
if (*minIndex<arr.Size()-1)
	*maxIndex = *minIndex+1;
else
	*maxIndex = *minIndex;
}


static double Interpol(double ene, float raeffG1, float raeffG2, float raeffE1, float raeffE2)
{
if (raeffE1==raeffE2)
	return raeffG1;
else if (raeffG1 == 0)
	return 0;
else if (raeffG2 == 0)
	return raeffG1 * exp(log(raeffE1/ene) / log(raeffE1));
else
	return exp(log(raeffG1) + log(raeffG2/raeffG1) * (log(ene/raeffE1) / log(raeffE2/raeffE1)));
}



/// Known issue: phi not used
double AeffGrid::Val(double ene, double theta, double phi) const
{
int m1, m2;
FindInterval(m_theta, theta, &m1, &m2);
int l1, l2;
FindInterval(m_energy, ene, &l1, &l2);

double a_inf_fit = Interpol(ene, m_aeffgrid(0,m1,l1), m_aeffgrid(0,m1,l2), m_energy[l1], m_energy[l2]);
double a_sup_fit = Interpol(ene, m_aeffgrid(0,m2,l1), m_aeffgrid(0,m2,l2), m_energy[l1], m_energy[l2]);

double aeff;
if (a_inf_fit==a_sup_fit)
	aeff = a_inf_fit;
else
	aeff = a_inf_fit  + ((a_sup_fit-a_inf_fit) / (m_theta[m2]-m_theta[m1])) * (theta - m_theta[m1]);
if (isnan(aeff))
	aeff = 0;

return aeff;
}





AeffGridAverage::AeffGridAverage(const char* aeffFileName)
 : AeffGrid(aeffFileName), m_hasEdp(false), m_emin(100), m_index(2.1)
{
m_emax = m_energy.Last();
m_avgValues.ResizeTo(m_aeffgrid.Dim(0), m_aeffgrid.Dim(1));
MakeGridAverage();
}

AeffGridAverage::AeffGridAverage(const char* aeffFileName, float emin, float emax, float index)
 : AeffGrid(aeffFileName), m_hasEdp(false), m_emin(emin), m_emax(emax), m_index(index)
{
m_avgValues.ResizeTo(m_aeffgrid.Dim(0), m_aeffgrid.Dim(1));
MakeGridAverage();
}

int AeffGridAverage::Read(const char* aeffFileName, float emin, float emax, float index)
{
int error = AeffGrid::Read(aeffFileName);
if (!error) {
	m_avgValues.ResizeTo(m_aeffgrid.Dim(0), m_aeffgrid.Dim(1));
	SetEIndex(emin, emax, index);
	}
return error;
}

int AeffGridAverage::LoadEdp(const char* edpFileName)
{
int error = 0;
if (edpFileName[0]=='N' && edpFileName[1]=='o' && edpFileName[2]=='n' && edpFileName[3]=='e' && edpFileName[4]==0)
	m_hasEdp = false;
else {
	error = m_edp.Read(edpFileName);
	m_hasEdp = !error;
	}
if (!error)
	MakeGridAverage();
return error;
}

/*
static float WeightedAvg(const VecD& eneArr, const float* aeffArr, double index)
{
if (index < 0)
	index = -index;
double result = aeffArr[0];
if (eneArr.Size()>1) {
	double specwt = 1.0;
	double coeff = 0;
	for (int i=0; i<eneArr.Size(); i++) {
		if (i+1==eneArr.Size())
			specwt *= pow(eneArr[i-1] / eneArr[i], 1.0-index);
		else
			specwt = pow(eneArr[i], 1.0-index) - pow(eneArr[i+1], 1.0-index) ;
		if (i==0)
			result = specwt * aeffArr[i];
		else
			result += specwt * aeffArr[i];
		coeff += specwt;
		}
	result /= coeff;
	}
return result;
}

int AeffGridAverage::MakeGridAverage()
{
int resultMask = 0;

int iMin = m_energy.GeomIndex(m_emin);
if (m_emin!=m_energy[iMin])
	resultMask = resultMask | 1;	/// Using different energy lower bound

int iMax = m_energy.Size()-1;
if (m_emax<=m_energy[iMax]) {
	iMax = m_energy.GeomIndex(m_emax);
	if (m_emax!=m_energy[iMax])
		resultMask = resultMask | 2;	/// Using different energy upper bound
	if (iMax>iMin)
		--iMax;
	}
else
	resultMask = resultMask | 4;	/// Upper bound treated as infinity

VecD evect(iMax-iMin+1);
for (int l=iMin ; l <= iMax; ++l)
	evect[l-iMin] = m_energy[l];

for (int i=0; i<m_avgValues.Dim(0); i++)
	for (int j=0; j<m_avgValues.Dim(1); j++)
		m_avgValues(i,j) = WeightedAvg(evect, &m_aeffgrid(i,j,iMin), m_index);
return resultMask;
}
*/



int AeffGridAverage::MakeGridAverage()
{
int resultMask = 0;

int iMin = m_energy.GeomIndex(m_emin);
if (m_emin!=m_energy[iMin])
	resultMask = resultMask | 1;	/// Using different energy lower bound

int eneChanCount = m_energy.Size();
int iMax = eneChanCount-1;
if (m_emax<=m_energy[iMax]) {
	iMax = m_energy.GeomIndex(m_emax);
	if (m_emax!=m_energy[iMax])
		resultMask = resultMask | 2;	/// Using different energy upper bound
	if (iMax>iMin)
		--iMax;
	}
else
	resultMask = resultMask | 4;	/// Upper bound treated as infinity

/// Calcolo del peso di ogni energia in base all'indice spettrale
VecF specwt(eneChanCount);
for (int i=0; i<eneChanCount-1; i++)
	specwt[i] = pow(double(m_energy[i]), 1.0-m_index) - pow(double(m_energy[i+1]), 1.0-m_index);
specwt[eneChanCount-1] = pow(double(m_energy[eneChanCount-1]), 1.0-m_index);

m_avgValues = 0.0f;
float normsum = 0;
for (int eobs = iMin; eobs <= iMax; eobs++)
	normsum += specwt[eobs];

int numtheta = m_avgValues.Dim(0);
int numphi = m_avgValues.Dim(1);
for (int thetaind = 0; thetaind < numtheta; thetaind++) {
	for (int phiind = 0; phiind < numphi; phiind++) {
		/// Calcolo della aeff da normalizzare
		VecF edpArr(eneChanCount);
		edpArr = 0.0f;
		float avgValue = 0.0f;
		for (int etrue = 0; etrue < eneChanCount; etrue++) {
			if (m_hasEdp)
				/// Calcolo della dispersione energetica totale per ogni canale di energia
				for (int eobs = iMin;  eobs <= iMax; eobs++)
					edpArr[etrue] += m_edp.Val(m_energy[etrue], m_energy[eobs], m_theta[thetaind], m_phi[phiind]);
			else
				edpArr[etrue] = (etrue<iMin || etrue>iMax) ? 0.0f : 1.0f;
			avgValue += edpArr[etrue] * specwt[etrue] * m_aeffgrid(thetaind, phiind, etrue);
			}
		m_avgValues(thetaind, phiind) = avgValue/normsum;
		}
	}
return resultMask;
}


/*
void AlikeAeffGridClass3::MakeWeightedAverageAeff(){

//	int eneChanCount = naxes[0];
//	int imin = eneChanCount - 1;
//	if ((*raeffenergy)[imin] >= emin) {
//	  for (imin=0; (*raeffenergy)[imin] < emin; imin++ );
//	  if (imin > 0 && emin / (*raeffenergy)[imin-1] < (*raeffenergy)[imin] / emin)
//		imin--;
//	}
//	int imax = eneChanCount - 1;
//	if ((*raeffenergy)[imax] >= emax) {
//	  for (imax=0; (*raeffenergy)[imax] < emax; imax++) ;
//	  if (imax > 0 && emax / (*raeffenergy)[imax-1] < (*raeffenergy)[imax] / emax)
//		imax--;
//	}
//	imin = imin;
//	imax = imax;
	const float* energies = *raeffenergy;
	int eneChanCount = naxes[0];
	naxesavg[0] = 1;
	naxesavg[1] = naxes[1];
	naxesavg[2] = naxes[2];
//	vector<Double_t> evect;
	imin = GeomIndex(energies, eneChanCount, emin);
	if (energies[imin]!=emin)
		cout << "WARNING: Lower bound for the energy channel of " << emin << " not found, using " << energies[imin] << " instead" << endl;
	imax = eneChanCount-1;
	if (emax<=energies[imax]) {
		imax = GeomIndex(energies, eneChanCount, emax);
		if (energies[imax]!=emax)
			cout << "WARNING: Upper bound for the energy channel of " << emax << " not found, using " << energies[imax] << " instead" << endl;
		if (imax>imin)
			--imax;
	}
	else
		cout << "Upper bound greater than the last channel, treated as infinity" << endl;

	int numphi = naxes[2];
	int numtheta = naxes[1];
	WeightedAverageAeff.ResizeTo(numphi, numtheta);
	WeightedAverageAeff = 0;

	if (index < 0)
		index = -index;
	
	/// Calcolo del peso di ogni energia in base all'indice spettrale
	TVectorD specwt(eneChanCount);
	for (int i=0; i < eneChanCount-1; i++)
	  specwt(i) =  pow((double) (*raeffenergy)[i] ,(double) (1.0 - index)) - pow((double)(*raeffenergy)[i+1], (double) (1.0 - index)) ;
	specwt(eneChanCount-1) = pow( (double) (*raeffenergy)[eneChanCount-1] , (double)(1.0 -index));
//	double normsum = 0;
//	for (int etrue = imin; etrue <= imax; etrue++)
//	  normsum += specwt(etrue);

	double normsum = 0;
	for (int eobs = imin;  eobs <= imax; eobs++)
		normsum += specwt(eobs);
		
	for (int thetaind = 0; thetaind < numtheta; thetaind++){
	  for (int phiind = 0; phiind < numphi; phiind++){
		  /// Calcolo della aeff da normalizzare
		  TVectorD edparr(eneChanCount);
		  for (int etrue = 0; etrue < eneChanCount; etrue++) {
			  /// Calcolo della dispersione energetica totale per ogni canale di energia
			for (int eobs = imin;  eobs <= imax; eobs++)
				edparr(etrue) += edpgrid.Val((*raeffenergy)[etrue], (*raeffenergy)[eobs], (*raefftheta)[thetaind], (*raeffphi)[phiind]);
			double energyDep = edparr(etrue) * specwt(etrue);
			WeightedAverageAeff(phiind, thetaind) += energyDep * raeffgrid[0][phiind][thetaind][etrue];
		  }
	    /// Normalizzazione dell'aeff
		  WeightedAverageAeff(phiind, thetaind) /= normsum;
	  }
	}
}

void AlikeAeffGridClass3::makegridavg(){

        imin = 0; 
	imax = 0;
	naxesavg[0] = 1;
	naxesavg[1] = naxes[1];
	naxesavg[2] = naxes[2];
	delete raeffgridavg;
	raeffgridavg = new  float*** [1];
	*raeffgridavg = new  float** [naxes[2]];
 	for (int i = 0; i < naxes[2]; i++) {
 		(*raeffgridavg)[i] = new float* [naxes[1]];
 		for (int j = 0; j < naxes[1]; j++) {
 			(*raeffgridavg)[i][j] = new float [1];	
 			}
 		}

	cout << "gridavg initialized" << endl;
	MakeWeightedAverageAeff();
	cout << "WeightedAverageAeff made" << endl;

 	for (int i = 0; i < naxes[2]; i++) {
 		for (int j = 0; j < naxes[1]; j++) {
 				(*raeffgridavg)[i][j][0] = WeightedAverageAeff(i,j);	
// 				cout << i << j << (*raeffgridavg)[i][j][0] << endl;
			}
		}
	cout << "gridavg made" << endl;
}

*/




































/** The initial version based on AlikeLib code
double AeffGridAverage::Valavg(double theta, double phi) const
{
int m1, m2;
FindInterval(m_theta, theta, &m1, &m2);

double a_inf_fit = Interpol(100, m_avgValues(0,m1), m_avgValues(0,m1), 100, 100);
double a_sup_fit = Interpol(100, m_avgValues(0,m2), m_avgValues(0,m2), 100, 100);

double aeff;
if (a_inf_fit==a_sup_fit)
	aeff = a_inf_fit;
else
	aeff = a_inf_fit  + ((a_sup_fit-a_inf_fit) / (m_theta[m2]-m_theta[m1])) * (theta - m_theta[m1]);
if (isnan(aeff))
	aeff = 0;
return aeff;
}
*/

/// Known issue: phi not used
double AeffGridAverage::AvgVal(double theta, double phi) const
{
int m1, m2;
FindInterval(m_theta, theta, &m1, &m2);
double inf = m_avgValues(0,m1);
double sup = m_avgValues(0,m2);
if (m_theta[m2]==m_theta[m1] || inf==sup)
	return inf;
return inf + ((sup-inf) / (m_theta[m2]-m_theta[m1])) * (theta-m_theta[m1]);
}


/// /////////////////////////////////////////////////////////////////////
///
/// class AeffEdpGridAverage
///
///
///
///

/**
int AeffEdpGridAverage::MakeGridAverage()
{
return AeffGridAverage::MakeGridAverage();
}
*/

/// /////////////////////////////////////////////////////////////////////
///
/// class PsfGrid
///
///
///
///

int PsfGrid::Read(const char* psfFileName)
{
FitsFile psfFile;
if (!psfFile.Open(psfFileName)) {
	cerr << "Error " << psfFile.Status() << " opening file" << psfFileName << endl;
	return psfFile.Status();
	}
psfFile.MoveAbsHDU(1);


psfFile.MoveAbsHDU(1);
VecL naxes;
if (GetImageSize(psfFile, naxes)) {
	m_psfgrid.ResizeTo(naxes[4], naxes[3], naxes[2], naxes[1], naxes[0]);
	long nelements = m_psfgrid.Size();
	int anynul = 0;
	fits_read_img_flt(psfFile, 0, 1, nelements, 0, m_psfgrid/*.Buffer()*/, &anynul, psfFile);
	psfFile.MoveAbsHDU(3);
	ReadFloatArray(psfFile, "RHO", m_psfrho);
	ReadFloatArray(psfFile, "PSI", m_psfpsi);
	ReadFloatArray(psfFile, "ENERGY", m_psfenergy);
	ReadFloatArray(psfFile, "POLAR_ANGLE", m_psftheta);
	ReadFloatArray(psfFile, "AZIMUTH_ANGLE", m_psfphi);
	}
return psfFile.Status();
}

/**
int PsfGrid::Read(const char* psfFileName)
{
int status = 0;
fitsfile* psfFile;
if (fits_open_file(&psfFile, psfFileName, READONLY, &status)) {
	cerr << "Error " << status << " opening file" << psfFileName << endl;
	return status;
	}

int naxis = 0;
long naxes[5];
/// int bitpix = 0; zzz unused
fits_movabs_hdu(psfFile, 1, NULL, &status);
/// fits_get_img_type(psfFile, &bitpix, &status); zzz unused
fits_get_img_dim(psfFile, &naxis, &status);
fits_get_img_size(psfFile, naxis, naxes, &status);

if (status==0) {

	m_psfgrid.ReshapeTo(5, naxes[4], naxes[3], naxes[2], naxes[1], naxes[0]);

	long nelements = m_psfgrid.Size();

	int anynul = 0;
	fits_read_img_flt(psfFile, 0, 1, nelements, 0, m_psfgrid.Buffer(), &anynul, &status);

	fits_movabs_hdu(psfFile, 3, NULL, &status);

	ReadFloatArray(psfFile, "RHO", &status, m_psfrho);
	ReadFloatArray(psfFile, "PSI", &status, m_psfpsi);
	ReadFloatArray(psfFile, "ENERGY", &status, m_psfenergy);
	ReadFloatArray(psfFile, "POLAR_ANGLE", &status, m_psftheta);
	ReadFloatArray(psfFile, "AZIMUTH_ANGLE", &status, m_psfphi);
	}

int closeStatus = 0;
fits_close_file(psfFile, &closeStatus);

return status;
}
*/

double PsfGrid::Val(double rho, double psi, double theta, double phi, double energy) const
{

int l = m_psfrho.LinearIndex(rho);
int m = m_psfpsi.LinearIndex(psi);
int n = m_psftheta.LinearIndex(theta);
int p = m_psfphi.LinearIndex(phi);
int q = m_psfenergy.GeomIndex(energy);
return m_psfgrid(p,n,q,m,l);
}







void PsfGridClass::Zero()
{
psfgrid = 0;
psfrho = psfpsi = psftheta = psfphi = psfenergy = 0;
naxes[0] = naxes[1] = naxes[2] =naxes[3] =naxes[4] = 0;
}



void PsfGridClass::Clean()
{
for (int i=0; i<naxes[4]; ++i) {
	for (int j=0; j<naxes[3]; ++j) {
		for (int k=0; k<naxes[2]; ++k) {
			for (int l=0; l<naxes[1]; ++l)
				delete[] psfgrid[i][j][k][l];
			delete[] psfgrid[i][j][k];
			}
		delete[] psfgrid[i][j];
		}
	delete[] psfgrid[i];
	}
delete[] psfgrid;
delete[] psfrho;
delete[] psfpsi;
delete[] psftheta;
delete[] psfphi;
delete[] psfenergy;
}



PsfGridClass::PsfGridClass(const char* psfFileName)
{
Zero();
Read(psfFileName);
}

bool PsfGridClass::Read(const char* psfFileName)
{
Clean();
int status = 0;
fitsfile * psfFile;
if (fits_open_file(&psfFile, psfFileName, READONLY, &status)) {
	cerr << "Error opening file" << psfFileName << endl;
	return status;
	}

int naxis = 0;
int bitpix = 0;
fits_movabs_hdu(psfFile, 1, NULL, &status);
fits_get_img_type(psfFile, &bitpix, &status);
fits_get_img_dim(psfFile, &naxis, &status);
fits_get_img_size(psfFile, naxis, naxes, &status);

float A[naxes[4]][naxes[3]][naxes[2]][naxes[1]][naxes[0]];

long nelements = naxes[0]*naxes[1]*naxes[2]*naxes[3]*naxes[4];

int anynul = 0;
fits_read_img_flt(psfFile, 0, 1, nelements, 0, A[0][0][0][0], &anynul, &status);

psfgrid = new float**** [naxes[4]];
for (int i = 0; i < naxes[4]; i++) {
	psfgrid[i] = new float*** [naxes[3]];
	for (int j = 0; j < naxes[3]; j++) {
		psfgrid[i][j] = new float** [naxes[2]];
		for (int k = 0; k < naxes[2]; k++) {
			psfgrid[i][j][k] = new float* [naxes[1]];
			for (int l = 0; l < naxes[1]; l++) {
				psfgrid[i][j][k][l] = new float [naxes[0]];
				}
			}
		}
	}

for (int i = 0; i < naxes[4]; i++) {
	for (int j = 0; j < naxes[3]; j++) {
		for (int k = 0; k < naxes[2]; k++) {
			for (int l = 0; l < naxes[1]; l++) {
				for (int n = 0; n < naxes[0]; n++) {
					psfgrid[i][j][k][l][n] = A[i][j][k][l][n];
					}
				}
			}
		}
	}

fits_movabs_hdu(psfFile, 3, NULL, &status);

long nrows = 0;
int colnum = 0;
int typecode = 0;
long width = 0;

char* colName = const_cast<char*>("RHO     ");
fits_get_colnum(psfFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(psfFile, colnum, &typecode, &nrows, &width, &status);
psfrho = new float[nrows];
fits_read_col(psfFile, TFLOAT, colnum, 1, 1, nrows, NULL, psfrho, NULL, &status);

colName = const_cast<char*>("PSI     ");
fits_get_colnum(psfFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(psfFile, colnum, &typecode, &nrows, &width, &status);
psfpsi= new float[nrows];
fits_read_col(psfFile, TFLOAT, colnum, 1, 1, nrows, NULL, psfpsi, NULL, &status);

colName = const_cast<char*>("ENERGY  ");
fits_get_colnum(psfFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(psfFile, colnum, &typecode, &nrows, &width, &status);
psfenergy = new  float [nrows];
fits_read_col(psfFile, TFLOAT, colnum, 1, 1, nrows, NULL, psfenergy, NULL, &status);

colName = const_cast<char*>("POLAR_ANGLE");
fits_get_colnum(psfFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(psfFile, colnum, &typecode, &nrows, &width, &status);
psftheta= new float[nrows];
fits_read_col(psfFile, TFLOAT, colnum, 1, 1, nrows, NULL, psftheta, NULL, &status);

colName = const_cast<char*>("AZIMUTH_ANGLE");
fits_get_colnum(psfFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(psfFile, colnum, &typecode, &nrows, &width, &status);
psfphi= new float[nrows];
fits_read_col(psfFile, TFLOAT, colnum, 1, 1, nrows, NULL, psfphi, NULL, &status);

int closeStatus = 0;
fits_close_file(psfFile, &closeStatus);

return status==0;
}


PsfGridClass::~PsfGridClass()
{
Clean();
}


double PsfGridClass::Val(double rho, double psi, double theta, double phi, double energy) const
{
int l = 0;
for (; rho>psfrho[l]; l++) ;
if (l>0 && rho-psfrho[l-1]<psfrho[l]-rho)
	--l;

int m = 0;
for (; psi>psfpsi[m]; ++m) ;
if (m>0 && psi-psfpsi[m-1]<psfpsi[m]-psi)
	--m;

int n = 0;
for (; theta>psftheta[n]; ++n) ;
if (n>0 && theta-psftheta[n-1]<psftheta[n]-theta)
	--n;

int p = 0;
for (; phi>psfphi[p]; ++p) ;
if (p>0 && phi-psfphi[p-1]<psfphi[p]-phi)
	--p;

int q = naxes[2]-1;
if (energy<psfenergy[naxes[2]-1]) {
	for (q=0; energy>psfenergy[q]; ++q) ;
	if (q > 0 && energy/psfenergy[q-1]<psfenergy[q]/energy)
		--q;
	}

return psfgrid[p][n][q][m][l];
}



#ifdef USE_TVECTOR
TVectorF PsfGridClass::energies() const
{
TVectorF ene(naxes[2]);
ene.SetElements(psfenergy);
return ene;
}



TVectorF PsfGridClass::rhos() const
{
TVectorF rho(naxes[0]);
rho.SetElements(psfrho);
return rho;
}

TVectorF PsfGridClass::psis() const
{
TVectorF psi(naxes[1]);
psi.SetElements(psfpsi);
return psi;
}
#endif




/// /////////////////////////////////////////////////////////////////////
///
/// class EdpGrid
///
///
///
///


double EdpGrid::Val(double trueE, double obsE, double theta, double phi) const
{
int l = m_edptrueenergy.GeomIndex(trueE);
int m = m_edpobsenergy.GeomIndex(obsE);
int n = m_edptheta.LinearIndex(theta);
int p = m_edpphi.LinearIndex(phi);
return m_edpgrid(p,n,m,l);
}




int EdpGrid::Read(const char* edpFileName)
{
FitsFile edpFile;
if (!edpFile.Open(edpFileName)) {
	cerr << "Error " << edpFile.Status() << " opening file" << edpFileName << endl;
	return edpFile.Status();
	}
edpFile.MoveAbsHDU(1);
VecL naxes;
if (GetImageSize(edpFile, naxes)) {
	m_edpgrid.ResizeTo(naxes[3], naxes[2], naxes[1], naxes[0]);
	long nelements = m_edpgrid.Size();
	int anynul = 0;
	fits_read_img_flt(edpFile, 0, 1, nelements, 0, m_edpgrid/*.Buffer()*/, &anynul, edpFile);
	edpFile.MoveAbsHDU(3);
	ReadFloatArray(edpFile, "TRUE_ENERGY", m_edptrueenergy);
	ReadFloatArray(edpFile, "OBS_ENERGY_CHANNEL", m_edpobsenergy);
	ReadFloatArray(edpFile, "POLAR_ANGLE", m_edptheta);
	ReadFloatArray(edpFile, "AZIMUTH_ANGLE", m_edpphi);
	}
return edpFile.Status();
}


/**
int EdpGrid::Read(const char* edpFileName)
{
int status = 0;
fitsfile* edpFile;
if (fits_open_file(&edpFile, edpFileName, READONLY, &status)) {
	cerr << "Error " << status << " opening file" << edpFileName << endl;
	return status;
	}

int naxis = 0;
long naxes[4];
/// int bitpix = 0; zzz unused
fits_movabs_hdu(edpFile, 1, NULL, &status);
/// fits_get_img_type(edpFile, &bitpix, &status); zzz unused
fits_get_img_dim(edpFile, &naxis, &status);
fits_get_img_size(edpFile, naxis, naxes, &status);

if (status==0) {

	m_edpgrid.ReshapeTo(4, naxes[3], naxes[2], naxes[1], naxes[0]);

	long nelements = m_edpgrid.Size();

	int anynul = 0;
	fits_read_img_flt(edpFile, 0, 1, nelements, 0, m_edpgrid.Buffer(), &anynul, &status);

	fits_movabs_hdu(edpFile, 3, NULL, &status);

	ReadFloatArray(edpFile, "TRUE_ENERGY", &status, m_edptrueenergy);
	ReadFloatArray(edpFile, "OBS_ENERGY_CHANNEL", &status, m_edpobsenergy);
	ReadFloatArray(edpFile, "POLAR_ANGLE", &status, m_edptheta);
	ReadFloatArray(edpFile, "AZIMUTH_ANGLE", &status, m_edpphi);
	}

int closeStatus = 0;
fits_close_file(edpFile, &closeStatus);

return status;
}
*/





void EdpGridClass::Zero()
{
edpgrid = 0;
edptrueenergy = edpobsenergy = edptheta = edpphi = 0;
naxes[0] = naxes[1] = naxes[2] = naxes[3] = 0;
}

void EdpGridClass::Clean()
{
for (int i=0; i<naxes[3]; ++i) {
	for (int j=0; j<naxes[2]; ++j) {
		for (int k=0; k<naxes[1]; ++k)
			delete[] edpgrid[i][j][k];
		delete[] edpgrid[i][j];
		}
	delete[] edpgrid[i];
	}
delete[] edpgrid;
delete[] edptrueenergy;
delete[] edpobsenergy;
delete[] edptheta;
delete[] edpphi;
}


bool EdpGridClass::Read(const char* edpFileName)
{
Clean();
int status = 0;
fitsfile* edpFile;
if (fits_open_file(&edpFile, edpFileName, READONLY, &status)) {
	cerr << "Error opening file" << edpFileName << endl;
	return status;
	}

int naxis = 0;
int bitpix = 0;
fits_movabs_hdu(edpFile, 1, NULL, &status);
fits_get_img_type(edpFile, &bitpix, &status);
fits_get_img_dim(edpFile, &naxis, &status);
fits_get_img_size(edpFile, naxis, naxes, &status);

float A[naxes[3]][naxes[2]][naxes[1]][naxes[0]];
long nelements = naxes[0]*naxes[1]*naxes[2]*naxes[3];

int anynul = 0;
fits_read_img_flt(edpFile, 0, 1, nelements, 0, A[0][0][0], &anynul, &status);

edpgrid = new float*** [naxes[3]];
for (int j = 0; j < naxes[3]; j++) {
	edpgrid[j] = new float** [naxes[2]];
	for (int k = 0; k < naxes[2]; k++) {
		edpgrid[j][k] = new float* [naxes[1]];
		for (int l = 0; l < naxes[1]; l++) {
			edpgrid[j][k][l] = new float [naxes[0]];
			}
		}
	}

for (int j = 0; j < naxes[3]; j++) {
	for (int k = 0; k < naxes[2]; k++) {
		for (int l = 0; l < naxes[1]; l++) {
			for (int n = 0; n < naxes[0]; n++) {
				edpgrid[j][k][l][n] = A[j][k][l][n];
				}
			}
		}
	}

fits_movabs_hdu(edpFile, 3, NULL, &status);
long nrows = 0;
int colnum = 0;
int typecode = 0;
long width = 0;
char* colName = const_cast<char*>("TRUE_ENERGY");
fits_get_colnum(edpFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(edpFile, colnum, &typecode, &nrows, &width, &status);
edptrueenergy = new float[nrows];
fits_read_col(edpFile, TFLOAT, colnum, 1, 1, nrows, NULL, edptrueenergy, NULL, &status);

colName = const_cast<char*>("OBS_ENERGY_CHANNEL");
fits_get_colnum(edpFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(edpFile, colnum, &typecode, &nrows, &width, &status);
edpobsenergy= new float[nrows];
fits_read_col(edpFile, TFLOAT, colnum, 1, 1, nrows, NULL, edpobsenergy, NULL, &status);

colName = const_cast<char*>("POLAR_ANGLE");
fits_get_colnum(edpFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(edpFile, colnum, &typecode, &nrows, &width, &status);
edptheta= new float[nrows];
fits_read_col(edpFile, TFLOAT, colnum, 1, 1, nrows, NULL, edptheta, NULL, &status);

colName = const_cast<char*>("AZIMUTH_ANGLE");
fits_get_colnum(edpFile, CASEINSEN, colName, &colnum, &status);
fits_get_eqcoltype(edpFile, colnum, &typecode, &nrows, &width, &status);
edpphi = new float[nrows];
fits_read_col(edpFile, TFLOAT, colnum, 1, 1, nrows, NULL, edpphi, NULL, &status);

int closeStatus = 0;
fits_close_file(edpFile, &closeStatus);

return status==0;
}


EdpGridClass::EdpGridClass(const char* edpFileName)
{
Zero();
Read(edpFileName);
}


EdpGridClass::~EdpGridClass()
{
Clean();
}


double EdpGridClass::Val(double E_true, double E_obs, double theta, double phi) const
{
/// return alikeEpd(edpgrid, edptrueenergy, edpobsenergy, edptheta, edpphi, E_true, E_obs, theta, phi, naxes);

/// cout << "trueenergy: " << trueenergy << " obsenergy: " << obsenergy << " theta: " << theta << " phi: " << phi << " naxes: " << naxes << endl;

int l = naxes[0]-1;
if (E_true < edptrueenergy[l]) {
	for ( l=0; E_true > edptrueenergy[l]; ++l) ;
	if (l>0 && (E_true / edptrueenergy[l-1] < edptrueenergy[l] / E_true))
		--l;
	}

int m = naxes[1]-1;
if (E_obs < edpobsenergy[naxes[1]-1]) {
	for (m=0 ; E_obs > edpobsenergy[m]; ++m) ;
	if (m > 0 && E_obs / edpobsenergy[m-1] < edpobsenergy[m] / E_obs)
		--m;
	}

int n = 0;
for (; theta>edptheta[n]; ++n) ;
if (n > 0 && theta - edptheta[n-1] < edptheta[n] - theta)
	--n;

int p = 0;
for (; phi>edpphi[p]; ++p) ;
if (p > 0 && phi - edpphi[p-1] < edpphi[p] - phi)
	--p;

/// zzz cerr << "old " << l << " " << m << " " << n << " " << p << endl;

return edpgrid[p][n][m][l];
}

#ifdef USE_TVECTOR
TVectorF EdpGridClass::trueenergies()
{
TVectorF ene(naxes[0]);
ene.SetElements(edptrueenergy);
return ene;
}

TVectorF EdpGridClass::obsenergies()
{
TVectorF ene(naxes[1]);
ene.SetElements(edpobsenergy);
return ene;
}
#endif


