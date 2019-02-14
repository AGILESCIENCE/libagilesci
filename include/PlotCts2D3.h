////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       Authors: Andrea Bulgarelli (INAF/IASF Bologna)
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

#ifndef _PLOT_CTS_2D_
#define _PLOT_CTS_2D_

void PlotCts2D(
	bool        shiftnorth_b,
	int         skysegmentation,
	const char* file,
	float       binsize,
	double      sourceradiusremove = 0,
	int         smoothing = 3,
	int         NCONREGSEARCH = 2,
	const char* outFile = "",
	bool        conregsearchType = 1,
	int         removeStep = 0,
	int         nstep = 1,
	int         bkgint = 20,
	const char* filesubtract = "",
	double      maxMapValue = -1);


void PlotCts2D_FINE(
	bool        shiftnorth_b,
	int         skysegmentation,
	const char* inputfile,
	float       binsize,
	double      sourceradiusremove = 0,
	int         smoothing = 3,
	int         NCONREGSEARCH = 2,
	const char* outFile = "",
	double      nearradious = 1,
	const char* expfile = "",
	double      minExp = 200,
	bool        conregsearchType = 1,
	int         removeStep = 0,
	int         nstep = 1,
	int         bkgint = 20,
	const char* filesubtract = "",
	double      maxMapValue = -1);




#endif

