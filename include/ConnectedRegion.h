////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       2008
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

#ifndef _CONNECTED_REGION_
#define _CONNECTED_REGION_


#include "TH2D.h"
#include "Labeling.h"

/**
@author Andrea Bulgarelli
*/

class ConnectedRegionList;

class ConnectedRegion /// : public TObject
{
	public:
		ConnectedRegion();

		~ConnectedRegion();

		void Print();

		void Reset();

		Int_t dim;

		Float_t center_x;

		Float_t center_y;

		Bool_t removed;

		Int_t first_bin_x; //indice M

		Int_t first_bin_y; //indice N

		Float_t n_b;

		Float_t sigma_x;

		Float_t sigma_y;

		//region
		Int_t M;
		Int_t N;
		Double_t* region;
		void SetRegion(Int_t k, Int_t M, Int_t N, Double_t *connected_regions);

		Int_t NumberOfContainedRegions(ConnectedRegionList* list);
		Int_t GetFirstContainedRegions(ConnectedRegionList* list);
		Bool_t Contains(ConnectedRegion* newregion);
		
		void CopyRegion(ConnectedRegion* region);

};

class ConnectedRegionList {
	public:
		ConnectedRegionList();

		~ConnectedRegionList();

		void RemoveNearPoints(Double_t radious);

		void Reset();

		void SetConnectedRegionsImage(Double_t* src);

		void SetOriginalImage(Double_t* src, TH2D* hist);

		void CalculateBaricenterAndError(Int_t nbindeg);

		void Print();

		void MergeRegions8Connection();
		/*3) ad ogni step mette in una lista blindata le regioni sicure. secondo le seguenti regole
		3.1) verifica se la regione X della lista blindata contiene 2 regioni della lista candidata. Se è così le rimuove entrambe dalla lista candidata
		3.2) di quelle rimaste dalla lista candidata, verifica per ogni blindata se è contenuta in una (ed una sola) candidata. In questo caso aggiorna la lista blindata con il candidato trovato (perchè la regione si è espansa) e la rimuove dalla lista candidata
		3.3) le superstiti della lista candidata sono tutte aggiunte alla blindata*/

		void MergeList(ConnectedRegionList* newlist);
	
		Int_t ncon;

		ConnectedRegion list[MAXLABELCONNECTED];

		TH2D* h2d_original;

		Double_t* connected_regions;

		Double_t* original_image;

		Int_t M;

		Int_t N;

		void SetIndexes();

		Int_t lastk;

		
};



#endif
