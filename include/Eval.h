////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       Authors: Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
//       Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
//       Tomaso Contessi (Nuove Idee sas)
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

#ifndef LIBAGILESCI_EVAL_H
#define LIBAGILESCI_EVAL_H

#include <vector>
#include "Intervals.h"
#include "AgileMap.h"

namespace eval
{

	struct MapspecEntry
	{
	    double fovradmin;
	    double fovradmax;
	    double emin;
	    double emax;
	    double index;
	};

	class Mapspec : public std::vector<MapspecEntry>
	{
	public:
	    Mapspec() {}
	    bool Store(const char* maplist);
	    void Print();
	};

	bool YTolTest(double ra_y, double dec_y, double ra_y0, double dec_y0, double y_tol);
	bool EarthTolTest(double earth_ra, double earth_dec, double earth_ra0, double earth_dec0, double earth_tol);
	bool AlbTest(double ra, double dec, double earth_ra, double earth_dec, double albrad);
	bool FovTest(Mapspec &maps, long i, double theta);

	bool LoadTimeList(const char* timelist, Intervals &intervals, double &tmin, double &tmax);

	int EvalExposure(const char *outfile, const char *sarFileName,
	                 const char *edpFileName, const char *maplist,
	                 const char *projection, double mdim, double mres,
	                 double la, double ba, double lonpole, double albrad,
	                 double y_tol, double roll_tol, double earth_tol,
	                 int phasecode, double binstep, int timestep,
	                 double index, double tmin, double tmax, double emin,
	                 double emax, double fovradmin, double fovradmax,
	                 const char *selectionFilename, const char *templateFilename,
	                 Intervals &intervals, std::vector< std::vector<double> > &exposures,
	                 bool saveMaps, bool sum_exposure, std::vector<double> &summed_exposures);

	int EvalCounts(const char *outfile, const char *projection, double tmin,
	               double tmax, double mdim, double mres, double la, double ba,
	               double lonpole, double emin, double emax, double fovradmax,
	               double fovradmin, double albrad, int phasecode, int filtercode,
	               const char *selectionFilename,  const char *templateFilename,
	               Intervals &intervals, std::vector< std::vector<int> > &counts,
	               bool saveMaps);

	/// Count the number of events (specified by radius) in a circle.
	/// The results are in counts
	int EvalCountsInRadius(const char *outfile, double tmin,
	               double tmax, double radius, double la, double ba,
	               double lonpole, double emin, double emax, double fovradmax,
	               double fovradmin, double albrad, int phasecode, int filtercode,
	               const char *selectionFilename,  const char *templateFilename,
	               Intervals &intervals, std::vector<int> &counts);

	int EvalGasMap(AgileMap &gasMap, AgileMap &expMap, const char* loresdiffuseFilename,
	               const char* hiresdiffuseFilename);

	int EvalGas(const char* outfile, const char* expfile, const char* diffusefile,
	            const char* hiresdiffusefile);

}
#endif
