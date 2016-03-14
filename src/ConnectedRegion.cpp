//
// C++ Implementation: ConnectedRegion
//
// Description: 
//
//
// Author: Andrea Bulgarelli <bulgarelli@iasfbo.inaf.it>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "ConnectedRegion.h"
#include "TMath.h"
#include "TH1D.h"
#include "TF1.h"
#include <iostream>
#include <cstdlib>


using namespace std;


ConnectedRegionList::ConnectedRegionList() {
	ncon = 0;
}

ConnectedRegionList::~ConnectedRegionList() {

}

void ConnectedRegionList::Reset() {
	for(int i=0; i<MAXLABELCONNECTED; i++) {
		list[i].Reset();
	}
}

void ConnectedRegionList::SetConnectedRegionsImage(Double_t* src) {
	connected_regions = new Double_t[M*N];
	for(int i=0; i< M; i++)
		for(int j=0; j< N; j++)
			connected_regions[i*M+j] = src[i*M+j];
}

void ConnectedRegionList::SetOriginalImage(Double_t* src, TH2D* hist) {
	original_image = new Double_t[M*N];
	for(int i=0; i< M; i++)
		for(int j=0; j< N; j++)
			original_image[i*M+j] = src[i*M+j];
	this->h2d_original = hist;
}

void ConnectedRegionList::MergeRegions8Connection() {
	SetIndexes();
	for(int k=0; k<lastk+1; k++) {
		if(list[k].dim == 1) {
			//cerco se nei quattro angoli 8 connessi ci sono regioni
			int i = list[k].center_x-1;
			int j = list[k].center_y-1;
// 			cout << "true center bin " << i << " " << j << " " << connected_regions[i * M + j] << endl;
// 			cout << connected_regions[(i+1) * M + j+1] << endl;
// 			cout << connected_regions[(i+1) * M + j-1] << endl;
// 			cout << connected_regions[(i-1) * M + j+1] << endl;
// 			cout << connected_regions[(i-1) * M + j-1] << endl;
			int newk = 0;
			if(connected_regions[(i+1) * M + j+1] != 0)
				connected_regions[i*M + j] = newk = connected_regions[(i+1) * M + j+1];
			if(connected_regions[(i+1) * M + j-1] != 0)
				connected_regions[i*M + j] = newk = connected_regions[(i+1) * M + j-1];
			if(connected_regions[(i-1) * M + j+1] != 0)
				connected_regions[i*M + j] = newk = connected_regions[(i-1) * M + j+1];
			if(connected_regions[(i-1) * M + j-1] != 0)
				connected_regions[i*M + j] = newk = connected_regions[(i-1) * M + j-1];
			if(newk != 0) {
				list[k].dim = 0;
				list[k].removed = true;
				list[newk-1].dim++;
				list[newk-1].region[i*M+j] = 1;
			}
// 			cout << newk << endl;
		}
	}
}

void ConnectedRegionList::RemoveNearPoints(Double_t radious) {
	cout << "#############################################" << endl;
	for(int k1=0; k1<MAXLABELCONNECTED; k1++) {
		if(list[k1].dim != 0)  
			for(int k2=k1+1; k2<MAXLABELCONNECTED; k2++) {
				if(list[k2].dim != 0) {
					if(TMath::Sqrt(TMath::Power(list[k1].center_x - list[k2].center_x, 2) + TMath::Power(list[k1].center_y - list[k2].center_y, 2)) <= radious) {
						cout << "REMOVE " << k2 << endl;
						list[k2].dim = 0;
					}
				}
			}
	}
}

void ConnectedRegionList::MergeList(ConnectedRegionList* newlist) {
	//3.1) verifica se la regione X della lista candidata contiene 2 regioni della lista blindata. Se è così le rimuove //la //regione della lista candidata
	SetIndexes();
	cout << "LASTK "<< lastk << " " << list << endl;
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(newlist->list[k].dim != 0)  
			if(newlist->list[k].NumberOfContainedRegions((ConnectedRegionList*)this) >= 2) {
				cout << "this new region " << k << " contains two old regions: removed " << endl;
 				newlist->list[k].dim = 0;
			}
	}
	//3.2) di quelle rimaste dalla lista candidata, verifica per ogni blindata se è contenuta in una (ed una sola) candidata. In questo caso aggiorna la lista blindata con il candidato trovato (perchè la regione si è espansa) e la rimuove dalla lista candidata
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(newlist->list[k].dim != 0)
			if(newlist->list[k].NumberOfContainedRegions((ConnectedRegionList*)this) == 1) {
				Int_t k1 = newlist->list[k].GetFirstContainedRegions((ConnectedRegionList*)this);
				cout << "region already exists in old list: update of old region " << k1 << endl;
				list[k1].CopyRegion( & (newlist->list[k]) );
				newlist->list[k].dim = 0;
			}
	}
	//3.3) le superstiti della lista candidata sono tutte aggiunte alla blindata
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(newlist->list[k].dim != 0) {
			cout << "copy " << k << endl;
			list[++lastk].CopyRegion(&newlist->list[k]);
			
		}
	}
	SetIndexes();
}

void ConnectedRegionList::CalculateBaricenterAndError(Int_t nbindeg) {
	Reset();
// 	SetIndexes();
	for(int k=1; k<=MAXLABELCONNECTED; k++)
		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++)
			{
				if(connected_regions[i*M + j] == k)
				{
					list[k-1].dim++;
					list[k-1].first_bin_x = i;
					list[k-1].first_bin_y = j;
				}
			}
	SetIndexes();
	//si setta la mappa della regione
	for(int k=1; k<=MAXLABELCONNECTED; k++) {
		if(list[k-1].dim != 0) {
			list[k-1].SetRegion(k, M, N, connected_regions);
		}
	}


	for(int k=0; k<MAXLABELCONNECTED; k++) list[k].center_x=list[k].center_y=list[k].n_b=0;
	for(int k=1; k<=MAXLABELCONNECTED; k++) {
		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++) 	{
				if(connected_regions[i*M + j] == k) {
					
					if(list[k-1].dim != 0) {
						
						list[k-1].n_b += original_image[i*M + j];
// 						cout << k << " " << (int) original_image[i*M + j] << " ";
						list[k-1].center_x += ( (i+1) * original_image[i*M + j] );
						list[k-1].center_y += ( (j+1) * original_image[i*M + j] );
					}
				}
			}
		
	}
// 	cout << endl;
	for(int k=0; k<MAXLABELCONNECTED; k++)
	{
		if(list[k].dim != 0) {
			list[k].center_x /= list[k].n_b;
			list[k].center_y /= list[k].n_b;
// 			cout << k << " " << list[k].center_x << " - " << list[k].center_y << "\n ";
		}
	}

		//calcolo dell'errore
	Double_t* src = new Double_t[M*N];
	TH2D *hist = (TH2D*)h2d_original->Clone("original2");
	Float_t radius_bin = nbindeg * 0.7; //prendo un raggio di N gradi attorno al peak della connected region
	for(int k=0; k<MAXLABELCONNECTED; k++)
	{
		if(list[k].dim != 0) {
			for(int i=0; i< M; i++)
				for(int j=0; j<N; j++)
					src[i*M + j] = 0;
			for(int i=1; i<= M; i++)
				for(int j=1; j<=N; j++)
				{
					Float_t dist = TMath::Sqrt( TMath::Power( list[k].center_x - i, 2 ) + TMath::Power( list[k].center_y - j, 2 ) );
					if(  dist  <=  radius_bin )
					{
						src[(i-1)*M + j-1] = h2d_original->GetBinContent(i,j);
					}
				}
	
			
			for(int i=0; i< M; i++)
				for(int j=0; j<N; j++)
					hist->SetBinContent(i+1, j+1, src[i*M + j]);
	
			TString name = "XCONREG_";
			name += (rand() * 1000000);
			TH1D* h1 = hist->ProjectionX(name, -1, -1, "d"); //mi metto in un raggio di 3 gradi
	
			//metodo 1 X:
			Float_t integralX = h1->Integral();
			Float_t sumX = h1->GetBinContent(list[k].center_x);
			Int_t nbinsX = 1;
			Int_t nstepX = 1;
			while(sumX < 	integralX * 0.683)
			{
				if(nbinsX % 2 == 0)
				{
					sumX += h1->GetBinContent(list[k].center_x + nstepX );
					nstepX++;
				}
				else
					sumX += h1->GetBinContent(list[k].center_x - nstepX );
				nbinsX++;
			}
			//esco che ho il numero di bins presi
			Float_t sigmaX_2 = ((nbinsX * integralX * 0.683) / sumX) / 2.0;
	
			//Metodo 2 X: fit con gaussiana
	
			h1->GetXaxis()->SetRangeUser(list[k].center_x - radius_bin * 1.1,  list[k].center_x + radius_bin * 1.1);
			h1->SetTitle(name);
			h1->Fit("gaus", "Q");
//			TF1* f1 = h1->GetFunction("gaus");
//			Float_t mean = f1->GetParameter(1); //mean gaus
//			Float_t sigma = f1->GetParameter(2); //sigma gaus
//			Float_t chi2 = f1->GetChisquare(); //chi square fit
//			Float_t ndf = f1->GetNDF(); //NDF fit
// 			cout << "X mean " << mean << " sigma " << sigma << " sigma in radius " << sigma / nbindeg << " ";
// 			cout << " X sigma 2 " << sigmaX_2 << " sigma2 in radius " << sigmaX_2 / nbindeg << endl;
			list[k].sigma_x = sigmaX_2 / nbindeg;
			//Y
	
			name = "YCONREG_";
			name += (rand() * 1000000);
			TH1D* h2 = hist->ProjectionY(name, -1, -1, "d"); //mi metto in un raggio di 3 gradi
	
			//metodo 1: Y
			Float_t integralY = h2->Integral();
			Float_t sumY = h2->GetBinContent(list[k].center_y);
			Int_t nbinsY = 1;
			Int_t nstepY = 1;
			while(sumY < 	integralY * 0.683)
			{
				if(nbinsY % 2 == 0)
				{
					sumY += h2->GetBinContent(list[k].center_y + nstepY );
					nstepY++;
				}
				else
					sumY += h2->GetBinContent(list[k].center_y - nstepY );
				nbinsY++;
			}
			//esco che ho il numero di bins presi
			Float_t sigmaY_2 = ((nbinsY * integralY * 0.683) / sumY ) / 2.0;
	
	
			//Metodo 2: fit con gaussiana Y
	
			h2->GetXaxis()->SetRangeUser(list[k].center_y - radius_bin * 1.1,  list[k].center_y + radius_bin * 1.1);
			h2->SetTitle(name);
			h2->Fit("gaus", "Q");
//			TF1* f2 = h2->GetFunction("gaus");
//			Float_t mean2 = f2->GetParameter(1); //mean gaus
//			Float_t sigma2 = f2->GetParameter(2); //sigma gaus
//			Float_t chi22 = f2->GetChisquare(); //chi square fit
//			Float_t ndf2 = f2->GetNDF(); //NDF fit
// 			cout << "Y mean " << mean2 << " sigma " << sigma2 << " sigma in radius " << sigma2 / nbindeg << " ";
// 			cout << " Y sigma 2 " << sigmaY_2 << " sigma2 in radius " << sigmaY_2 / nbindeg << endl;
			list[k].sigma_y = sigmaY_2 / nbindeg;
		}
	}
	delete []src;
}

void ConnectedRegionList::SetIndexes() {
	lastk = 0;
	ncon = 0;
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(list[k].dim != 0) {
			lastk = k;
			ncon++;
		}
	}

}


void ConnectedRegionList::Print() {
	SetIndexes();
	for(int k=0; k<lastk+1; k++) {
		cout << k << " ";
		list[k].Print();
	}
}

ConnectedRegion::ConnectedRegion()
{
	removed = false;
	dim = 0;
	region = 0;
}


ConnectedRegion::~ConnectedRegion()
{
}

void ConnectedRegion::SetRegion(Int_t k, Int_t M, Int_t N, Double_t *connected_regions) {
	this->M = M;
	this->N = N;
	if(region == 0)
		region = new Double_t[M*N];
	else {
		for(int i=0; i<M; i++)
			for(int j=0; j<N; j++)
				region[i*M+j] = 0;
	}
	for(int i=0; i<M; i++)
		for(int j=0; j<N; j++)
			if(connected_regions[i*M+j] == k)
				region[i*M+j] = 1;
}

void ConnectedRegion::Print() {
	cout << dim << " " << center_x << " +/- " << sigma_x << " " << center_y << " +/- " << sigma_y << " "  << first_bin_x << " " << first_bin_y  << " bin: ";
	for(int i=0; i<M; i++)
		for(int j=0; j<N; j++)
			if(region[i*M+j] != 0)
				cout << "[" << i << ", " << j << "]";
	cout << endl;
}

void ConnectedRegion::Reset() {
	dim = center_x = center_y = removed = first_bin_x = first_bin_y = n_b = sigma_x = sigma_y = 0;
}

Int_t ConnectedRegion::NumberOfContainedRegions(ConnectedRegionList* oldlist) {
// 	cout << "pointer of oldlist " << oldlist << endl;
	Int_t num = 0;
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(oldlist->list[k].dim != 0) {
// 			cout << "compare "; Print();
//  			cout << "with "; oldlist->list[k].Print();
			if(Contains(&oldlist->list[k]))
				num++;
		}
	}
// 	cout << "NUM "<< num << endl;
	return num;
}

Int_t ConnectedRegion::GetFirstContainedRegions(ConnectedRegionList* oldlist) {
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(oldlist->list[k].dim != 0) {
			if(Contains(&oldlist->list[k]))
				return k;
		}
	}
	return -1;
}



Bool_t ConnectedRegion::Contains(ConnectedRegion* oldregion) {
	Bool_t contains = true;
	for(int i=0; i<M; i++)
		for(int j=0; j<N; j++) {
			if(oldregion->region[i*M+j] != 0) {
				if(region[i*M+j] == 0) {
					contains = false;
					break;
				}
					
			}
		}
	return contains;
}


void ConnectedRegion::CopyRegion(ConnectedRegion* region) {
		

		/// if (dim == 0 || gm.AlgorithmType == 1 ) { zzz gm.AlgorithmType was not initialized
		if (dim == 0) {
			this->center_x = region->center_x;

			this->center_y = region->center_y;
		}

		this->dim = region->dim;

		this->removed= region->removed;

		this->first_bin_x= region->first_bin_x; //indice M

		this->first_bin_y= region->first_bin_y; //indice N

		this->n_b= region->n_b;

		this->sigma_x= region->sigma_x;

		this->sigma_y= region->sigma_y;

		//region
		this->M= region->M;
		this->N= region->N;
		if(this->region == 0)
			this->region = new Double_t[M*N];
		if(region->region != 0)
			for(int i=0; i<M ; i++)
				for(int j=0; j<N; j++)
					this->region[i*M+j] = region->region[i*M+j];
}
