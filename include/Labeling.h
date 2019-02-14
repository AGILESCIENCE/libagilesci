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

#ifndef _LABELING_
#define _LABELING_



const int MAXLABELCONNECTED = 2000;
const int BG = 0;

void SetImage(int dim_x, int dim_y, double* space_density_1);
/// void SetSpaceDensity2D(Int_t i, Int_t j, Double_t value);
/// Double_t GetSpaceDensity2D(Int_t i, Int_t j);
/// Int_t SecondStep(unsigned int last_i, unsigned int NewLabel, int* C );
int LabelingDuePassi();
int LabelingDuePassi8();	/// Unused for now

#endif
