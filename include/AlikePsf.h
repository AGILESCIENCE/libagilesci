/*
 * Copyright (c) 2005-2017 
 *     AGILE Team
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 * GPL License
*/


#ifndef _ALIKE_PSF_
#define _ALIKE_PSF_

#include <string>

#include "CalibUtils.h"
#include "AgileMap.h"


class AlikeCircle
{
public:
	AlikeCircle(): m_map() {}
	AlikeCircle(double lCenter, double bCenter, double rAnal, const AgileMap& map)
	 : m_map() { Reset(lCenter, bCenter, rAnal, map); }
	AlikeCircle(const AlikeCircle& another): m_map(another.m_map) {}
	~AlikeCircle() {}

	AlikeCircle& operator=(const AlikeCircle &another) { m_map=another.m_map; return *this; }

	void Reset(double lCenter, double bCenter, double rAnal, const AgileMap& map);

	int GetNoElements() const { return m_map.Size(); }
	int GetBinNo(int i) const { return m_map[i]; }
	int operator[](int i) const { return m_map[i]; }

private:
	VecI m_map;
};

/**
class AlikeCircle
{
public:
	AlikeCircle(): m_map(0), m_count(0) {}
	AlikeCircle(double lCenter, double bCenter, double rAnal, const AgileMap& map)
	 : m_map(0) { Reset(lCenter, bCenter, rAnal, map); }
	AlikeCircle(const AlikeCircle& another) : m_map(0) { Copy(another); }
	~AlikeCircle() { delete[] m_map; }

	AlikeCircle& operator=(const AlikeCircle &another) { return Copy(another); }

	void Reset(double lCenter, double bCenter, double rAnal, const AgileMap& map);

	int GetNoElements() const { return m_count; }
	int GetBinNo(int i) const { return m_map[i]; }
	int operator[](int i) const { return m_map[i]; }

private:
	AlikeCircle& Copy(const AlikeCircle &another);

private:
	int*   m_map;
	int    m_count;
};
*/


class AlikeNorm
{
public:
	AlikeNorm(): m_eInf(100), m_eSup(50000), m_normFactor(1), m_integratortype(1) {}
	AlikeNorm(const AlikeNorm& o): m_eInf(o.m_eInf), m_eSup(o.m_eSup), m_normFactor(o.m_normFactor), m_integratortype(o.m_integratortype) {}
	~AlikeNorm() {}

	AlikeNorm& operator=(const AlikeNorm& o) { m_eInf=o.m_eInf, m_eSup=o.m_eSup, m_normFactor=o.m_normFactor; m_integratortype=o.m_integratortype; return *this; }

	double GetEnergyInf() const { return m_eInf; }
	double GetEnergySup() const { return m_eSup; }
	double GetNormFactor() const { return m_normFactor; }
	double PLnuFnu(double eMin, double eMax, double index, int type = 0);
	double PLExpCutOffnuFnu(double eMin, double eMax, double index, double m_par2, int type = 0);
	double LogParabolanuFnu(double eMin, double eMax, double index, double m_par2, double m_par3, int type = 0);
	double PLSuperExpCutOffnuFnu(double eMin, double eMax, double index, double m_par2, double m_par3, int type = 0);
	
protected:
	/// WARNING: After calling SetEnergyRange the user should call UpdateNorm()
	void SetEnergyRange(double eInf, double eSup) { m_eInf = eInf; m_eSup = eSup; }
	double UpdateNorm(double eMin, double eMax, double index, bool norm=1);
	double UpdateNormPLExpCutOff(double eMin, double eMax, double index, double par2, bool norm=1);
	double UpdateNormLogParabola(double eMin, double eMax, double index, double par2, double par3, bool norm=1);
	double UpdateNormPLSuperExpCutOff(double eMin, double eMax, double index, double par2, double par3, bool norm=1);
	
	
	int m_integratortype;
	
	double UpdateIntegrator(TF1& f, double eMin, double eMax, double eInf, double eSup, bool norm);
	
private:
	double m_eInf;
	double m_eSup;
	double m_normFactor;
	
};




class AlikePsfTables: private AeffGrid, private PsfGrid, private EdpGrid
{
public:
	AlikePsfTables(): m_sinRhoArr(1), m_deltaRho(1), m_psfFileName(""), m_raeffFileName(""), m_edpFileName("") {}
	~AlikePsfTables() {}

public:
	bool Read(const char* psfFileName, const char* raeffFileName, const char* edpFileName);

	const char* PsfFileName() const { return m_psfFileName.c_str(); }
	const char* RaeffFileName() const { return m_raeffFileName.c_str(); }
	const char* EdpFileName() const { return m_edpFileName.c_str(); }

	double AeffVal(double ene, double theta, double phi) const
		{ return AeffGrid::Val(ene, theta, phi); }
	double PsfVal(double rho, double psi, double theta, double phi, double energy) const
		{ return PsfGrid::Val(rho, psi, theta, phi, energy); }
	double EdpVal(double E_true, double E_obs, double theta, double phi) const
		{ return EdpGrid::Val(E_true, E_obs, theta, phi); }

	int PsfeCount() const { return PsfGrid::m_psfenergy.Size(); }
	float PsfEnergy(int i) const { return PsfGrid::m_psfenergy[i]; }

	float DeltaRho() const { return m_deltaRho; }
	int   RhoCount() const { return PsfGrid::m_psfrho.Size(); }
	float Rho(int i) const { return PsfGrid::m_psfrho[i]; }
	float SinRho(int i) const { return m_sinRhoArr[i]; }

	const VecF& PsfEnerges() const { return PsfGrid::m_psfenergy; }
	const VecF& Rhos() const { return PsfGrid::m_psfrho; }
	const VecF& SinRhos() const { return m_sinRhoArr; }
	
	const VecF& EdpTrueEnerges() const { return EdpGrid::m_edptrueenergy; }
	const VecF& EdpObsEnerges() const { return EdpGrid::m_edpobsenergy; }
	const VecF& EdpTheta() const { return EdpGrid::m_edptheta; }
	const VecF& EdpPhi() const { return EdpGrid::m_edpphi; }
	const Mat4F& EdpVal() const { return EdpGrid::m_edpgrid; }

private:	/// Data
	VecF        m_sinRhoArr;
	float       m_deltaRho;
	std::string m_psfFileName;
	std::string m_raeffFileName;
	std::string m_edpFileName;
};




class AlikePsfSource: public AgileMap, public AlikeNorm
{
public:	/// Construction
	AlikePsfSource(): AgileMap(), AlikeNorm(), m_psfTab(0), m_theta(0), m_srcL(0), m_srcB(0), m_index(2.1), m_typefun(0), m_par2(3000.0), m_par3(2.0), m_psfArr(0), m_edpArr(0), m_specwt(0), m_init(false), m_init_index(2.1), m_init_par2(3000.0), m_init_par3(2.0) {}
	~AlikePsfSource() {}
	void Set(
		const AlikePsfTables* psfTab,
		const AgileMap &inmap,
		double eInf, double eSup,
		double theta,
		double srcL, double srcB, double index, int typefun, double par2, double par3);

public:	/// Getting info
	const AlikePsfTables* PsfTab() const { return m_psfTab; }

	double GetSrcL() const { return m_srcL; }
	double GetSrcB() const { return m_srcB; }
	double GetIndex() const { return m_index; }
	double GetTypeFun() const { return m_typefun; }
	double GetPar2() const { return m_par2; }
	double GetPar3() const { return m_par3; }

	double DegDistance(double l, double b) const { return SphDistDeg(m_srcL, m_srcB, l, b); }

	const double* PsfArr() const { return m_psfArr; }
	int GetPsfCount() const { return m_psfTab->RhoCount(); }

public:	/// Operations
	enum Changes { NoChanges=0, IndexChanged=1, PositionChanged=2, IndexPositionChanged=3 };
	Changes SetSrcData(double srcL, double srcB, double index, double par2, double par3, bool force=false);

protected:	/// Data
	const AlikePsfTables* m_psfTab;
	double m_theta;

	double m_srcL;
	double m_srcB;
	double m_index;
	int m_typefun;
	double m_par2;
	double m_par3;
	double m_init_index;
	double m_init_par2;
	double m_init_par3;
	bool m_init;

	VecD   m_psfArr;

private:	/// Ancillary data
	VecD   m_edpArr;
	VecD   m_specwt;

/// Not allowed
private:
	AlikePsfSource(const AlikePsfSource& another);
	AlikePsfSource& operator=(const AlikePsfSource& another);
	void Copy(const AlikePsfSource& another);
};




#endif
