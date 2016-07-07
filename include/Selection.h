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

#ifndef LIBAGILESCI_SELECTION_H
#define LIBAGILESCI_SELECTION_H

#include <string>
#include <Intervals.h>

namespace selection
{

std::string TimesExprString(const Intervals& intvs);

std::string LogExprString(const Intervals& intvs, int phasecode, int timeStep);

int MakeSelection(const char *fileList, Intervals& selection,
                  std::string expr, const char *selectionFilename,
                  const char *templateFilename);

}

#endif
