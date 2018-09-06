/*
 * Copyright (c) 2005-2016
 *     Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *     Tomaso Contessi (Nuove Idee sas)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 * GPL License
*/

#ifndef LIBAGILESCI_EVAL_H
#define LIBAGILESCI_EVAL_H

#include <vector>
#include "Intervals.h"
#include "AgileMap.h"

namespace eval
{
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
                 bool saveMaps);

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
