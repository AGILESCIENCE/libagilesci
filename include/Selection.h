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

#ifndef LIBAGILESCI_SELECTION_H
#define LIBAGILESCI_SELECTION_H

#include <string>
#include <Intervals.h>

namespace selection
{

std::string TimesExprString(const Intervals& intvs);

std::string LogExprString(const Intervals& intvs, int phasecode, int timeStep);

std::string EvtExprString(const Intervals &intvs, double emin, double emax,
                          double albrad, double fovradmax, double fovradmin,
                          int phasecode, int filtercode);

int MakeSelection(const char *fileList, Intervals& selection,
                  std::string expr, const char *selectionFilename,
                  const char *templateFilename);

}

#endif
