/*
 * Copyright (c) 2005-2016
 *     Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *     Tomaso Contessi (Nuove Idee sas)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

#ifndef LIBAGILESCI_EVAL_H
#define LIBAGILESCI_EVAL_H

#include <vector>
#include <Intervals.h>

namespace eval
{

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

}

#endif
