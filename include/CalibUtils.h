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
//       14/Feb/2012
//             First release: V1.0
//             Author: Andrew Chen, Tomaso Contessi (IASF-Milano), Andrea Bulgarelli (IASF-Bologna)
//             Based on 3 classes previously defined in AlikeLib.h
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       All rights reserved.
////////////////////////////////////////////////////////////////////////////////////

#ifndef _CALIB_UTILS_
#define _CALIB_UTILS_


#include "ArrayOf.h"


typedef ArrayOf<4, float>  Mat4F;
typedef ArrayOf<5, float>  Mat5F;



class PsfGrid
{
public:
	PsfGrid() {}
	PsfGrid(const char* psfFileName) { Read(psfFileName); }
	
	int Read(const char* psfFileName); /// Return the CFITSIO error
	double Val(double rho, double psi, double theta, double phi, double energy) const;

	const Mat5F& Values() const { return m_psfgrid; }
	const VecF& Energies() const { return m_psfenergy; }
	const VecF& Rhos() const { return m_psfrho; }
	const VecF& Psis() const { return m_psfpsi; }
	const VecF& Thetas() const { return m_psftheta; }
	const VecF& Phis() const { return m_psfphi; }

protected:
	Mat5F m_psfgrid;
	VecF  m_psfphi; 		/// dim(0)
	VecF  m_psftheta; 	/// dim(1)
	VecF  m_psfenergy;	/// dim(2)
	VecF  m_psfpsi; 		/// dim(3)
	VecF  m_psfrho; 		/// dim(4)
};




class AeffGrid
{
public:
	AeffGrid() {}
	AeffGrid(const char* aeffFileName) { Read(aeffFileName); }
	~AeffGrid() {}

	int Read(const char* aeffFileName); /// Return the CFITSIO error

	double Val(double ene, double theta, double phi) const;
	
	const  Mat3F& Values() const { return m_aeffgrid; }
	const  VecF& Theta() const { return m_theta; }
	const  VecF& Phi() const { return m_phi; }
	const  VecF& Energies() const { return m_energy; }

protected:
	Mat3F m_aeffgrid;
	VecF  m_theta;
	VecF  m_phi;
	VecF  m_energy;
};


class EdpGrid
{
public:
	EdpGrid() {}
	EdpGrid(const char* edpFileName) { Read(edpFileName); }

	int Read(const char* edpFileName);

	double Val(double trueE, double obsE, double theta, double phi) const;

	const Mat4F& Values() const { return m_edpgrid; }

	const VecF& TrueEnergies() const { return m_edptrueenergy; }
	const VecF& ObsEnergies() const { return m_edpobsenergy; }
	const VecF& Thetas() const { return m_edptheta; }
	const VecF& Phis() const { return m_edpphi; }

protected:
	Mat4F m_edpgrid;
	VecF  m_edptrueenergy;
	VecF  m_edpobsenergy;
	VecF  m_edptheta;
	VecF  m_edpphi;
};




class AeffGridAverage: public AeffGrid
{
public:
	AeffGridAverage(): AeffGrid(), m_edp(), m_hasEdp(false), m_emin(100), m_emax(1000), m_index(2.1), m_avgValues() {}
	AeffGridAverage(const char* aeffFileName);
	AeffGridAverage(const char* aeffFileName, float emin, float emax, float index);

	/// Return the CFITSIO error
	int Read(const char* aeffFileName, float emin, float emax, float index);
	int LoadEdp(const char* edpFileName);

	const MatF& AvgValues() const { return m_avgValues; }

	double AvgVal(double theta, double phi) const;

	double Emin() const { return m_emin; }
	double Emax() const { return m_emax; }
	double Index() const { return m_index; }

	void SetEmin(double emin) { m_emin = emin; MakeGridAverage(); }
	void SetEmax(double emax) { m_emax = emax; MakeGridAverage(); }
	void SetE(double emin, double emax) { m_emin = emin; m_emax = emax; MakeGridAverage(); }
	void SetIndex(double index) { m_index = index>0?index:-index; MakeGridAverage(); }
	void SetEIndex(double emin, double emax, double index) { m_emin = emin; m_emax = emax; m_index = index>0?index:-index; MakeGridAverage(); }

protected:
	EdpGrid m_edp;
	bool    m_hasEdp;
	float   m_emin;
	float   m_emax;
	float   m_index;
	MatF    m_avgValues;
	int MakeGridAverage();
};


/**
class AeffEdpGridAverage: public AeffGridAverage, private EdpGrid
{
public:
	AeffEdpGridAverage(const char* aeffFileName, const char* edpFileName)
	 : AeffGridAverage(aeffFileName), EdpGrid(edpFileName) {}
	AeffEdpGridAverage(const char* aeffFileName, const char* edpFileName, float emin, float emax, float index)
	 : AeffGridAverage(aeffFileName, emin, emax, index), EdpGrid(edpFileName) {}

protected:
	virtual int MakeGridAverage();
};
*/



/// Old classes to delete


class PsfGridClass
{
public:
	PsfGridClass() { Zero(); }
	PsfGridClass(const char* psfFileName);
	~PsfGridClass();
	
	bool Read(const char* psfFileName);
	int getAxes(int n) const { return naxes[n]; }
	double Val(double rho, double psi, double theta, double phi, double energy) const;

	float ValOf(int i1, int i2, int i3, int i4, int i5) const { return psfgrid[i1][i2][i3][i4][i5]; }
	int DimOf(int axis) const { return naxes[axis]; }

#ifdef USE_TVECTOR
	TVectorF energies() const;
	TVectorF rhos() const;
	TVectorF psis() const;
#else
	int PsfRhoCount() const { return naxes[0]; }
	const float* PsfRhoArr() const { return psfrho; }
	int PsfEnergyCount() const { return naxes[2]; }
	const float* PsfEnergyArr() const { return psfenergy; }
	/// int PsfPsiCount() const { return naxes[1]; }
	/// const float* PsfPsiArr() const { return psfpsi; }
#endif

private:
	void Zero();
	void Clean();
private:
	long       naxes[5];
	float***** psfgrid;
	float*     psfrho;		/// naxes[0]
	float*     psfpsi;		/// naxes[1]
	float*     psftheta;		/// naxes[3]
	float*     psfphi;		/// naxes[4]
	float*     psfenergy;	/// naxes[2]
};






class EdpGridClass
{
public:
	EdpGridClass() { Zero(); }
	EdpGridClass(const char* edpFileName);
	~EdpGridClass();


	float ValOf(int i1, int i2, int i3, int i4) const { return edpgrid[i1][i2][i3][i4]; }
	int DimOf(int axis) const { return naxes[axis]; }

	bool Read(const char* edpFileName);
	int getAxes(int n) const { return naxes[n]; }
	double Val(double E_true, double E_obs, double theta, double phi) const;
#ifdef USE_TVECTOR
	TVectorF trueenergies();
	TVectorF obsenergies();
#else
	int EdpTrueEnergyCount() const { return naxes[0]; }
	const float* EdpTrueEnergyArr() const { return edptrueenergy; }
	int EdpObsEnergyCount() const { return naxes[1]; }
	const float* EdpObsEnergyArr() const { return edpobsenergy; }
#endif
private:
	void Zero();
	void Clean();
private:
	long      naxes[4];
	float**** edpgrid;
	float*    edptrueenergy;	/// naxes[0]
	float*    edpobsenergy;		/// naxes[1]
	float*    edptheta;
	float*    edpphi;
};



#endif
