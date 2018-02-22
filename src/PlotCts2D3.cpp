


#include "PlotCts2D3.h"



#define DRAW 1

#include "TNamed.h"
#include "TVirtualFitter.h"
#include "TH1.h"
#include "TF1.h"
/// #include "TMinuit.h"

#include "TH2.h"
#include "TCanvas.h"
// #include "TSpectrum.h"
#include "TMarker.h"

// #include "TMatrixD.h"




#include "Labeling.h"
#include "ConnectedRegion.h"
#include "TMath.h"

/// #include "AlikeMap4.h"
#include "AgileMap.h"
#define AlikeMap AgileMap
#define GetNrows Rows
#define GetNcols Cols
#define AlikeSphdistDeg SphDistDeg
#include <iostream>
#include <fstream>
using namespace std;



#define DIMRES 9
#define NPEKAMAXFOREACHCUT 20



static Double_t x_b[MAXLABELCONNECTED];
static Double_t y_b[MAXLABELCONNECTED];
static Double_t n_b[MAXLABELCONNECTED];
static Double_t sigma_x[MAXLABELCONNECTED];
static Double_t sigma_y[MAXLABELCONNECTED];
static Double_t SNR[MAXLABELCONNECTED];
static Double_t pixel[MAXLABELCONNECTED];


static Int_t nconnectedregion_found = 0;



void ResetArrays() {
	for(int i=0; i<MAXLABELCONNECTED; i++) {
		x_b[i] = y_b[i] = n_b[i] = sigma_x[i] = sigma_y[i] = SNR[i] = pixel[i] = 0;
	}
	nconnectedregion_found = 0;
}

void gaussian(float* kernel, int r)
{
	int rr = 2*r+1;
	int ksz = rr*rr;
	float sigma = r/2.;

	float kt = 0;
	float a = 1./(sigma*sigma);
	float c = 1./(sigma*sigma);
	for (int y=-r; y<=r; y++)
	{
		for (int x=-r; x<=r; x++)
		{
			if ((x*x + y*y) <= r*r)
			{
				float v = exp(-.5*(a*x*x + c*y*y));
				kernel[(y+r)*rr+(x+r)] = v;
				kt += v;
			}
		}
	}

	// normalize kernel
	for (int aa=0; aa<ksz; aa++)
		kernel[aa] /= kt;
}

void convolve(float* kernel, float* src, float* dest,
              int width, int height, int r)
{
	int rr = 2*r+1;

	float* dptr = dest;
	for (int jj=0; jj<height; jj++)
	{
		for (int ii=0; ii<width; ii++, dptr++)
		{

			for (int nn=jj-r, qq=0; nn<=jj+r; nn++, qq++)
			{
				if (nn>=0 && nn<height)
				{
					int nd = nn*width;
					int qd = qq*rr;
					for (int mm=ii-r, pp=0; mm<=ii+r; mm++, pp++)
					{
						if (mm>=0 && mm<width)
							*dptr += src[nd+mm]*kernel[qd+pp];
					}
				}
			}
		}
	}
}

void RemoveNoise(Int_t nstep, TH2D* hist, Int_t M, Int_t N, Int_t bkgint)
{
	//rimozione del background
	Int_t nsteps = nstep;
	Int_t nbins = M / nsteps;
	cout << "REMOVE NOISE: nsteps " << nsteps << " - nbins " << nbins << endl;
	for(int i=0; i<nsteps; i++)
	{
		Int_t base = nbins*i+1;
		TString name = "X";
		name += i;
		Int_t y1, y2;
		if(i < nsteps - 1)
		{
			y1 = base;
			y2 = base+nbins-1;
		}
		else
		{
			y1 = base;
			y2 = M;
		}
		cout << "bin step " << y1 << " " << y2 << endl;
		TH1D* h1 = hist->ProjectionX(name, y1, y2);
		h1->SetTitle(name);

		new TCanvas; h1->Draw();

		/// TH1D* h1b =  (TH1D*)h1->ShowBackground(bkgint);		compiler says its unused
		Int_t div = y2-y1+1;
		cout << "bin of background: " << h1->GetNbinsX() <<". / " << div << endl;
		for(int y=y1; y<=y2; y++)
			for(int x=1; x<=N; x++)
			{
				Double_t diff = hist->GetBinContent(x,y) - h1->GetBinContent(x) / (div);
				diff = diff > 0 ? diff : 0;
				hist->SetBinContent(x, y, diff);
			}

	}
}

void Smooth(TH2D* hist, Int_t M, Int_t N, int r1)
{
	int r=r1;
	int rr = 2*r+1;
	int ksz = rr*rr;

	float* kernel = new float[ksz];

	float* src = new float[M*N];
	float* dest = new float[M*N];
	for(int i=0; i< M; i++)
		for(int j=0; j< N; j++)
			src[i*M+j] = dest[i*M+j]=0;

	for(int i=1; i<= M; i++)
		for(int j=1; j<=N; j++)
			src[(i-1)*M + j-1] = hist->GetBinContent(i,j);
	//copy
	// 	for(int i=0; i< M; i++)
	// 		for(int j=0; j<N; j++)
	// 			dest[i*M + j] = src[i*M+j];

	gaussian(kernel, r);
	convolve(kernel, src, dest, M, N, r);
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			hist->SetBinContent(i+1, j+1, dest[i*M + j]); //M-1 per effettuare la rotazione lungo l'asse verticale
	delete[] kernel;
	delete[] src;
	delete[] dest;
}



/**   UNUSED FUNCTION

Double_t* SearchProfileY(TH2D* hist, Int_t M, Int_t N, Int_t smoothing, Int_t& index)
{
	index = 0;
	Int_t nn = N-smoothing+1;
	cout << nn << endl;
	Double_t* res = (Double_t*) new Double_t[nn*DIMRES * NPEKAMAXFOREACHCUT];

	for(int k=1; k<=nn; k++)
	{
		TString name = "A";
		name += k;
		TH1D* h1 = hist->ProjectionY(name, k, k+smoothing-1);


		// 		cout << k << ", " << k+smoothing-1 << " - " << index*DIMRES + 0 << " " << index << endl;
		TSpectrum* sp = new TSpectrum;
		Int_t npeaks = sp->Search(h1);
		Float_t* peaksX = sp->GetPositionX();
		Float_t* peaksY = sp->GetPositionY();
		for(int i=0; i<npeaks; i++)
		{
			h1->GetXaxis()->SetRangeUser(peaksX[i]-smoothing, peaksX[i]+smoothing);
			h1->Fit("gaus", "Q");
			TF1* f1 = h1->GetFunction("gaus");
			//  			cout << "PEAKX: " << peaksX[i] << " Y: " << peaksY[i] << " - mean " << f1->GetParameter(1) << " +/- " << f1->GetParameter(2) << " - chi2 " << f1->GetChisquare() << "/" << f1->GetNDF() << endl;

			res[index*DIMRES + 0] = k;
			res[index*DIMRES + 1] = k+smoothing-1;
			res[index*DIMRES + 2] = peaksX[i]; //punto del picco su asse X
			res[index*DIMRES + 3] = peaksY[i]; //valore dell'histo nel picco
			res[index*DIMRES + 4] = f1->GetParameter(1); //mean gaus
			res[index*DIMRES + 5] = f1->GetParameter(2); //sigma gaus
			res[index*DIMRES + 6] = f1->GetChisquare(); //chi square fit
			res[index*DIMRES + 7] = f1->GetNDF(); //NDF fit
			res[index*DIMRES + 8] = i; //number of peak for this cut
			index++;

		}

		delete sp;

	}
	return res;
}

*/


/**   UNUSED FUNCTION

Double_t* SearchProfileX(TH2D* hist, Int_t M, Int_t N, Int_t smoothing, Int_t& index)
{
	index = 0;
	Int_t nn = M-smoothing+1;
	cout << nn << endl;
	Double_t* resX = (Double_t*) new Double_t[nn*DIMRES * NPEKAMAXFOREACHCUT];

	for(int k=1; k<=nn; k++)
	{
		TString name = "B";
		name += k;
		TH1D* h1 = hist->ProjectionX(name, k, k+smoothing-1);
		// 		cout << k << ", " << k+smoothing-1 << " - " << index*DIMRES + 0 << " " << index << endl;
		TSpectrum* sp2 = new TSpectrum;
		Int_t npeaks = sp2->Search(h1);
		Float_t* peaksX = sp2->GetPositionX();
		Float_t* peaksY = sp2->GetPositionY();
		for(int i=0; i<npeaks; i++)
		{
			h1->GetXaxis()->SetRangeUser(peaksX[i]-smoothing, peaksX[i]+smoothing);
			h1->Fit("gaus", "Q");
			TF1* f1 = h1->GetFunction("gaus");
			//  			cout << "PEAKX: " << peaksX[i] << " Y: " << peaksY[i] << " - mean " << f1->GetParameter(1) << " +/- " << f1->GetParameter(2) << " - chi2 " << f1->GetChisquare() << "/" << f1->GetNDF() << endl;

			resX[index*DIMRES + 0] = k;
			resX[index*DIMRES + 1] = k+smoothing-1;
			resX[index*DIMRES + 2] = peaksX[i]; //punto del picco su asse X
			resX[index*DIMRES + 3] = peaksY[i]; //valore dell'histo nel picco
			resX[index*DIMRES + 4] = f1->GetParameter(1); //mean gaus
			resX[index*DIMRES + 5] = f1->GetParameter(2); //sigma gaus
			resX[index*DIMRES + 6] = f1->GetChisquare(); //chi square fit
			resX[index*DIMRES + 7] = f1->GetNDF(); //NDF fit
			resX[index*DIMRES + 8] = i; //number of peak for this cut
			index++;

		}

		delete sp2;

	}
	return resX;
}
*/

/**    UNUSED FUNCTION

void DrawImage(TString file)
{
	new TCanvas();
	TFShortImg * img2 = (TFShortImg *)TFReadImage(file, NULL, 1);
	cout << "img2 " << img2 << endl;
	UInt_t size[2];
	img2->GetSize(size);
	Short_t* da = img2->GetDataArray();
	cout << da << endl;
	for (UInt_t x = 0; x < size[1]; x++)
		for (UInt_t y = 0; y < size[0]; y++)
		{
			(*img2)[y][x] = 1.0;//hist->GetBinContent(x+1, y+1);
			// 			if(da[x*size[1] + y])
			cout << x << " " << y << " " << da[x*size[1] + y] << endl;
		}
	// create the image display for the just read and slightly
	// modified image
	AstroImage * disp = new AstroImage(img2);
	// draw the image
	disp->Draw();
	disp->StartPaletteEditor();
}
*/


//la seguente funzione cerca crea una immagine con un numero di regioni connesse pari a quelle passate come argomento, a partire
//dall'histo passato come argomento
//typeconreg = true significa che mi fermo quando trovo nconreg (prendo quindi solo la punta della nuova conreg)
//typeconreg = false mi fermo quando trovo nconreg + 1
//nbindeg = numero di bin che sono contenuti in un grado

TH2D* CreateConnectedRegions(TH2D* hist, Int_t nconreg, Int_t typeconreg, Int_t nbindeg, Int_t skysegmentation)
{
	hist->Scale(100/hist->GetMaximum()); //scale da 0 a 100
	Int_t M = hist->GetNbinsX(); //dim X
	Int_t N = hist->GetNbinsY();
	
	Int_t max = hist->GetMaximum() + 1; //si arrotonda
	cout << "--------------------------------" << endl;
	cout << "max " << max << endl;

	Double_t* src = new Double_t[M*N];
	Double_t* src_tmp = new Double_t[M*N];
	Double_t* src_conreg_final = new Double_t[M*N];
	for(int i=0; i< M; i++)
		for(int j=0; j< N; j++)
			src[i*M+j] = 0;

	// 	for(int i=1; i<= M; i++)
	// 		for(int j=1; j<=N; j++) {
	// 			src[(i-1)*M + j-1] = hist->GetBinContent(i,j);
	// 	}

	SetImage(M, N, src);
	Int_t last_nconreg = 0;
	int ncon = 0; //number of connected regions
	for(int k=max; k>=0; k--)
	{
// 		cout << "K: " << k << endl;
		for(int i=1; i<= M; i++)
			for(int j=1; j<=N; j++)
			{
				if( ((Int_t)hist->GetBinContent(i,j)) == k)
				{
					src[(i-1)*M + j-1] = hist->GetBinContent(i,j);
					// 					cout << "F: " << i << "," << j  << " " << src[(i-1)*M + j-1] << endl;
				}
			}
		//per ogni k, provo a calcolare il numero di connessioni connesse

		ncon = LabelingDuePassi();

		cout << " [" << k << ", " << ncon << "] ";
		//-----------------------------------------------------------------------------------------------------------------
		//metodo 1: per cercare regioni connesse. Mi fermo quando ne ho trovate in numero pari a quanto specificato in input
		//-----------------------------------------------------------------------------------------------------------------
		if(typeconreg == 1) //mi fermo appena trovo n con reg o piu'
			if(ncon >= nconreg)
				break;

		//-----------------------------------------------------------------------------------------------------------------
		//metodo 2: cerca fino a quando non trovo un numero di regione connesse maggiore di 1 di quelle specificate,
		//oppure mi fermo quando il numero di regioni connesse inizia a diminuire (pesco background)
		//-----------------------------------------------------------------------------------------------------------------
		if(typeconreg == 0)
		{ 	//mi fermo appena trovo n+1 con reg e prendo n con reg
			if(ncon > nconreg || last_nconreg > ncon)
			{ 	//se ho un numero maggiore di regioni connesse di quelle richieste, oppure
				//se il numero di regioni connesse diminuisce, allora mi fermo
				for(int i=0; i< M; i++)
					for(int j=0; j<N; j++)
					{
						src[i*M + j] = src_tmp[i*M + j];
					}
				ncon = LabelingDuePassi();
				break;
			}
			else
				for(int i=0; i< M; i++)
					for(int j=0; j<N; j++)
					{
						src_tmp[i*M + j] = src[i*M + j];
					}
		}
		last_nconreg = ncon;

	}
	cout << endl;
	//*******************
	//**draw result
	TH2D *bori3 = (TH2D*)hist->Clone("end step 1");
	TString histname = "ES1";
	histname += skysegmentation;
	bori3->SetName(histname);

	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src[i * M + j] > 0)
				bori3->SetBinContent(i+1, j+1, src[i * M + j]); //hist->GetBinContent(i+1, j+1)
			else
				bori3->SetBinContent(i+1, j+1, 0);
	bori3->SetTitle("con reg found - number");
	#ifdef DRAW
	new TCanvas("con reg found number", "con reg found number");
	bori3->Draw("LEGO2");
	#endif
	TH2D *bori4 = (TH2D*)hist->Clone("end step 1");
	histname = "ESS2";
	histname += skysegmentation;
	bori4->SetName(histname);

	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src[i * M + j] > 0)
				bori4->SetBinContent(i+1, j+1, hist->GetBinContent(i+1, j+1));
			else
				bori4->SetBinContent(i+1, j+1, 0);
	bori4->SetTitle("con reg found");
	#ifdef DRAW
	new TCanvas("con reg found", "con reg found");
	bori4->Draw("LEGO2");
	#endif

	//****************************************************************
	//correzione delle regioni connesse trovate in base al rapporto SNR
	//****************************************************************
	//1) calcola il numero di fotoni di background per ogni pixel
	//1.1) sottrai le regioni connesse trovate
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src[i*M + j] == 0)
				src_tmp[i * M + j] = hist->GetBinContent(i+1, j+1);
			else
				src_tmp[i * M + j] = 0;
	//1.2) e ora calcola il valor medio
	Double_t mean_bkg = 0, points_bkg = 0;
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src_tmp[i * M + j] != 0)
			{
				points_bkg++;
				mean_bkg += src_tmp[i * M + j];
			}
	mean_bkg /= points_bkg;
	cout << "mean value of the background " << mean_bkg << endl;

	//*******************
	//**draw result
	TH2D *bori2 = (TH2D*)hist->Clone("end step 2");
	histname = "ES2";
	histname += skysegmentation;
	bori2->SetName(histname);

	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			bori2->SetBinContent(i+1, j+1, src_tmp[i * M + j]);

	bori2->SetTitle("con reg subtracted");
	#ifdef DRAW
	new TCanvas("con reg subtracted", "con reg subtracted");
	bori2->Draw("COL");
	#endif

	//****************************************************************
	//estrazione di una regione connessa e calcolo di tutto!
	//****************************************************************
	cout << "Number of connected regions found: " << ncon << endl;
	// 	Int_t nreg = 1;
	for(int nreg = 1; nreg<=ncon; nreg++)
	{
		Double_t ncon_sig = 0, ncon_pixel = 0;
		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++)
				if(src[i*M + j] == nreg)
				{
					src_tmp[i * M + j] = hist->GetBinContent(i+1, j+1);
					ncon_pixel++;
					ncon_sig += src_tmp[i * M + j];
				}
				else
					src_tmp[i * M + j] = 0;

		//si calcola il numero di pixel e l'integrale
		Double_t noise = ncon_pixel * mean_bkg;
		Double_t SNR = ncon_sig / TMath::Sqrt(ncon_sig + 2 * noise);
		cout << "************ connected region " << nreg << endl;
		cout << "SNR = " << SNR << endl;
		cout << "pixel = " << ncon_pixel << endl;
		cout << "integral of noise = " << noise << endl;
		cout << "integral of signal = " << ncon_sig << endl;
		
	}
	//****************************************************************
	//********** CALCOLO DEL BARICENTRO - METHOD 1: Calculate center of mass della regione connessa
	//****************************************************************

	//calcola il baricentro di ogni regione connessa
	for(int k=0; k<ncon; k++) x_b[k]=y_b[k]=n_b[k]=0;
	for(int k=1; k<=ncon; k++)
		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++)
			{
				if(src[i*M + j] == k)
				{
					n_b[k-1] += hist->GetBinContent(i+1, j+1);
					x_b[k-1] += ((i+1) * hist->GetBinContent(i+1, j+1));
					y_b[k-1] += ((j+1) * hist->GetBinContent(i+1, j+1));
				}
			}
	cout << "baricenter calculation......" << endl;

	//salva immagine originale
	TH2D *boriginal = (TH2D*)hist->Clone("original");
	histname = "ORI1";
	histname += skysegmentation;
	boriginal->SetName(histname);

	//salva istogramma derivante dall'estrazione delle regioni connesse
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src[i*M + j] == 0)
				hist->SetBinContent(i+1, j+1, 0);
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			hist->SetBinContent(i+1, j+1, src[i*M + j]);
	TH2D *bconreg = (TH2D*)hist->Clone("connected_regions");
	histname = "CRR1";
	histname += skysegmentation;
	bconreg->SetName(histname);

	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			src_conreg_final[i*M + j] = 0;

	//aggiungo anche la regione attorno al peak della connected region
	Float_t radius_bin = nbindeg * 0.7; //prendo un raggio di N gradi attorno al peak della connected region
	cout << "radius_bin " << radius_bin << endl;
	#ifdef DRAW
	new TCanvas;
	#endif
	for(int k=0; k<ncon; k++)
	{
		x_b[k] /= n_b[k];
		y_b[k] /= n_b[k];
		cout << x_b[k] << " - " << y_b[k] << " ";

		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++)
				src[i*M + j] = 0;

		for(int i=1; i<= M; i++)
			for(int j=1; j<=N; j++)
			{
				Float_t dist = TMath::Sqrt( TMath::Power( x_b[k] - i, 2 ) + TMath::Power( y_b[k] - j, 2 ) );
				if(  dist  <=  radius_bin )
				{
					src[(i-1)*M + j-1] = boriginal->GetBinContent(i,j);
					src_conreg_final[(i-1)*M + j-1] = boriginal->GetBinContent(i,j);
				}
			}


		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++)
				hist->SetBinContent(i+1, j+1, src[i*M + j]);

		TString name = "XCONREG_";
		name += k;
		TH1D* h1 = hist->ProjectionX(name); //mi metto in un raggio di 3 gradi


		//metodo 1 X:
		Float_t integralX = h1->Integral();
		Float_t sumX = h1->GetBinContent(x_b[k]);
		Int_t nbinsX = 1;
		Int_t nstepX = 1;
		while(sumX < 	integralX * 0.683)
		{
			if(nbinsX % 2 == 0)
			{
				sumX += h1->GetBinContent(x_b[k] + nstepX );
				nstepX++;
			}
			else
				sumX += h1->GetBinContent(x_b[k] - nstepX );
			nbinsX++;
		}
		//esco che ho il numero di bins presi
		Float_t sigmaX_2 = ((nbinsX * integralX * 0.683) / sumX) / 2.0;

		//Metodo 2 X: fit con gaussiana

		h1->GetXaxis()->SetRangeUser(x_b[k] - radius_bin * 1.1,  x_b[k] + radius_bin * 1.1);
		h1->SetTitle(name);
		// 		h1->Draw();
		h1->Fit("gaus", "Q");
		TF1* f1 = h1->GetFunction("gaus");
		Float_t mean = f1->GetParameter(1); //mean gaus
		Float_t sigma = f1->GetParameter(2); //sigma gaus
		Float_t chi2 = f1->GetChisquare(); //chi square fit
		Float_t ndf = f1->GetNDF(); //NDF fit
		cout << "X mean " << mean << " sigma " << sigma << " sigma in radius " << sigma / nbindeg << " ";
		cout << " X sigma 2 " << sigmaX_2 << " sigma2 in radius " << sigmaX_2 / nbindeg << endl;
		sigma_x[k] = sigmaX_2 / nbindeg;
		//Y

		name = "YCONREG_";
		name += k;
		TH1D* h2 = hist->ProjectionY(name); //mi metto in un raggio di 3 gradi

		//metodo 1: Y
		Float_t integralY = h2->Integral();
		Float_t sumY = h2->GetBinContent(y_b[k]);
		Int_t nbinsY = 1;
		Int_t nstepY = 1;
		while(sumY < 	integralY * 0.683)
		{
			if(nbinsY % 2 == 0)
			{
				sumY += h2->GetBinContent(y_b[k] + nstepY );
				nstepY++;
			}
			else
				sumY += h2->GetBinContent(y_b[k] - nstepY );
			nbinsY++;
		}
		//esco che ho il numero di bins presi
		Float_t sigmaY_2 = ((nbinsY * integralY * 0.683) / sumY ) / 2.0;

		//Metodo 2: fit con gaussiana Y

		h2->GetXaxis()->SetRangeUser(y_b[k] - radius_bin * 1.1,  y_b[k] + radius_bin * 1.1);
		h2->SetTitle(name);
		// 		new TCanvas;
		// 		h2->Draw();
		h2->Fit("gaus", "Q");
		TF1* f2 = h2->GetFunction("gaus");
		Float_t mean2 = f2->GetParameter(1); //mean gaus
		Float_t sigma2 = f2->GetParameter(2); //sigma gaus
		Float_t chi22 = f2->GetChisquare(); //chi square fit
		Float_t ndf2 = f2->GetNDF(); //NDF fit
		cout << "Y mean " << mean2 << " sigma " << sigma2 << " sigma in radius " << sigma2 / nbindeg << " ";
		cout << " Y sigma 2 " << sigmaY_2 << " sigma2 in radius " << sigmaY_2 / nbindeg << endl;
		sigma_y[k] = sigmaY_2 / nbindeg;
	}


	// 	hist->Draw("COL");
	bconreg->SetTitle("Connected region");
	#ifdef DRAW
	new TCanvas;
	bconreg->Draw("COL");
	#endif
	for(int k=0; k<ncon; k++)
	{
		TMarker* m = new TMarker(x_b[k], y_b[k], 2);
		m->Draw();
	}


	
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			hist->SetBinContent(i+1, j+1, src_conreg_final[i*M + j]);
	hist->SetTitle("Extended connected region");
	#ifdef DRAW
	new TCanvas;
	hist->Draw("COL");
	#endif

	nconnectedregion_found = ncon;

    delete []src_tmp;
    delete []src_conreg_final;
	return hist;

}

bool SetSkyRegion(
	TH2D* hist,
	int skysegmentation,
	const AlikeMap& map,	/// TFFloatImg* img,
	double sourceradiusremove,
	double lcenter,
	double bcenter)
{
if(skysegmentation == 0)
	return true;

/// AstroImage2 * disp = new AstroImage2(img);
int M = hist->GetNbinsX();
int N = hist->GetNbinsY();
bool hasValues = false;
for(int i=0; i< M; i++)
	for(int j=0; j<N; j++) {
		double l = map.l(i, j);
		double b = map.b(i, j);
/**
		Double_t b = disp->Get_b(i+1, j+1);	
		Double_t l = disp->Get_l(i+1, j+1);
		Double_t dist = disp->Sphdistance(lcenter, bcenter, l, b, 1);
*/
		if (skysegmentation == 1 && b < 10)
			hist->SetBinContent(i+1, j+1, 0);
		else if (skysegmentation == 2 && (b >= 10 || b <= -10))
			hist->SetBinContent(i+1, j+1, 0);
		else if (skysegmentation == 3 && b > -10)
			hist->SetBinContent(i+1, j+1, 0);
		else {
			double dist = AlikeSphdistDeg(l, b, lcenter, bcenter);
			if (dist > sourceradiusremove + 1.0)
				hist->SetBinContent(i+1, j+1, 0);
			else
				hasValues = true;
			}
		}
/// delete disp;
return hasValues;
}


/*skysegmentation
	0 -> all the sky
	1 -> b > 10
	2 ->  -10 <= b <= 10
	3 -> b < -10
*/


static TH2D MakeTH2D(const AlikeMap& map, const char* name="", const char* title="")
{
int nRows = map.GetNrows();
int nCols = map.GetNcols();
TH2D h(name, title, nCols, 0.5, 0.5+nCols, nRows, 0.5, 0.5+nRows);
h.SetDirectory(0);
for (int i=0; i<nRows; ++i)
	for (int j=0; j<nCols; ++j)
		h.SetCellContent(i+1, j+1, map(i,j));
return h;
}



extern void PlotCts2D(
	bool        shiftnorth_b,
	int         skysegmentation,
	const char* file,
	float       binsize,
	double      sourceradiusremove,
	int         smoothing,
	int         NCONREGSEARCH,
	const char* outFile,
	bool        conregsearchType,
	int         removeStep,
	int         nstep,
	int         bkgint,
	const char* filesubtract,
	double      maxMapValue)
{
	ResetArrays();
	/// gStyle->SetPalette(1);	zzz ‘gStyle’ was not declared in this scope
	Int_t nbindeg = 1/binsize;
	cout << "nbindeg " << nbindeg << endl;

	AlikeMap map(file);
   double lcenter = map.GetMapCenterL();
   double bcenter = map.GetMapCenterB();

	TH2D hist = MakeTH2D(map);
	TH2D hist2(hist);
	if(SetSkyRegion(&hist, skysegmentation, map, sourceradiusremove, lcenter, bcenter) == 0)
		return;
	if(SetSkyRegion(&hist2, skysegmentation, map, sourceradiusremove, lcenter, bcenter) == 0)
		return;

	int M = hist.GetNbinsX(); //dim X
	int N = hist.GetNbinsY();
	//l'histo 2D come numerazione bin va da (1,1) a (M,N)
	cout << "M*N = " << M << " * " << N << endl;

	if (filesubtract[0]) {
		AlikeMap mapSub(filesubtract);
      /// TH2D histSub(mapSub);
		TH2D histSub = MakeTH2D(mapSub);
		if (SetSkyRegion(&histSub, skysegmentation, mapSub, sourceradiusremove, lcenter, bcenter) == false)
			return;
		for (int i=0; i< M; i++)
			for (int j=0; j<N; j++) {
				Float_t sub = TMath::Abs(hist.GetBinContent(i+1, j+1) - hist2.GetBinContent(i+1, j+1));
				sub = sub > 0 ? sub : 0;
				hist.SetBinContent(i+1, j+1, sub);
				}
		}

	if (maxMapValue != -1)
		for (int i=0; i< M; i++)
			for (int j=0; j<N; j++)
				if (hist.GetBinContent(i+1, j+1) > maxMapValue) {
					cout << "A " << hist.GetBinContent(i+1, j+1) << endl;
					hist.SetBinContent(i+1, j+1, 0);
					hist2.SetBinContent(i+1, j+1, 0);
					}

	hist2.SetTitle("Original");
	#ifdef DRAW
 	new TCanvas("Original", "Original");
 	hist2.Draw("COL"); //NON CANCELLARE!!!!
	#endif

	if (smoothing)
		Smooth(&hist, M, N, smoothing);

	TH2D bori(hist);
	bori.SetNameTitle("original", "Smoothed");
	#ifdef DRAW
	new TCanvas("Smoothed", "Smoothed");
	bori.Draw("COL");
	#endif

	for (int i=0; i< removeStep; i++) {
		cout << "RemoveNoise" << endl;
		RemoveNoise(nstep, &hist, M, N, bkgint);
		}
// 	Int_t dimY;
// 	Int_t dimX;
	//    	Double_t* resY = SearchProfileY(hist, M, N, smoothing, dimY);
	//  	Double_t* resX = SearchProfileX(hist, M, N, smoothing, dimX);

	/*
		for(int i=0; i<dimY; i++) {
	  		cout << "Ymin " << resY[i*DIMRES+0] << " Ymax " << resY[i*DIMRES+1] << " Npeak " << resY[i*DIMRES+8] << " PEAK_X: " << resY[i*DIMRES+2] << " Y_value: " << resY[i*DIMRES+3] << " - mean " << resY[i*DIMRES+4] << " +/- " << resY[i*DIMRES+5] << " - chi2 " << resY[i*DIMRES+6] << "/" << resY[i*DIMRES+7] << endl;
		}
		for(int i=0; i<dimX; i++) {
	  		cout << "Xmin " << resX[i*DIMRES+0] << " Xmax " << resX[i*DIMRES+1] << " Npeak " << resX[i*DIMRES+8] << " PEAK_Y: " << resX[i*DIMRES+2] << " X_value: " << resX[i*DIMRES+3] << " - mean " << resX[i*DIMRES+4] << " +/- " << resX[i*DIMRES+5] << " - chi2 " << resX[i*DIMRES+6] << "/" << resX[i*DIMRES+7] << endl;
		}
	*/

/// zzz Tom: hh2 is returned but not used
	TH2D* hh2 = CreateConnectedRegions(&hist, NCONREGSEARCH, conregsearchType, nbindeg, skysegmentation);  //0.5 deg per ogni bin - 2 bin contengono un grado
// 	new TCanvas();
// 	hh2->Draw("COL");

	// 	return hist;
	//provo a salvare l'histo in un img
// 	new TCanvas();
// 	TFShortImg * img2 = (TFShortImg *)TFReadImage(file, NULL, 1);
// 	cout << "img2 " << img2 << endl;
// 	UInt_t size[2];
// 	img2->GetSize(size);
// 	for (UInt_t x = 0; x < size[1]; x++)
// 		for (UInt_t y = 0; y < size[0]; y++)
// 		{
			//  			(*img2)[y][x] = hist->GetBinContent(x+1, y+1);
			// 			cout << y << " " << x << endl;
// 		}
	// create the image display for the just read and slightly
	// modified image
 	new TCanvas;
// 	AstroImage2 * disp = new AstroImage2(img);
	// draw the image
 	/// disp->Draw();		zzz only statement where disp is used (but now its not supported)
	// 	disp->StartPaletteEditor();
	std::ofstream asciiFile1;
	std::ofstream asciiFile2;

	TString listFile(outFile);
	listFile += ".list";
	TString regFile(outFile);
	regFile += ".reg";


	if(skysegmentation == 0 || skysegmentation == 1) {
		asciiFile1.open(listFile);
		asciiFile2.open(regFile);
		
	} else {
		asciiFile1.open(listFile, std::ios::app);
		asciiFile2.open(regFile, std::ios::app);
	}
	asciiFile2 << "global color=green font=\"helvetica 9 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\ngalactic\n";
// 	std::ofstream asciiFile1(outFile + ".list");
// 	std::ofstream asciiFile2(outFile + ".reg");
// 	asciiFile2 << "global color=green font=\"helvetica 10 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\ngalactic\n";
	
// 	Double_t lcenter = disp->Get_l(bincenter_x, bincenter_y);
// 	Double_t bcenter = disp->Get_b(bincenter_x, bincenter_y) - (shiftnorth_b?binsize:0.0);

	cout << "CENTER OF THE MAP " << lcenter << " " << bcenter << endl;	

	static int kold=0;
	int knew = 0;
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		knew = k;
// 		cout << x_b[i] << " " << y_b[i] << " ";
		Double_t lcoo = map.l(x_b[k]-1, y_b[k]-1) ;
		Double_t bcoo = map.b(x_b[k]-1, y_b[k]-1) ;

		cout << lcoo << " " << bcoo << endl;
		Double_t dista = AlikeSphdistDeg(lcenter, bcenter, lcoo, bcoo);
		if(sourceradiusremove == 0 || dista <= sourceradiusremove) {
			
			asciiFile1 << "1 0+/-0 1 1 1 " << map.l(x_b[k], y_b[k]) << " " << map.b(x_b[k], y_b[k]) - (shiftnorth_b?binsize:0.0) << " 1 1 0.00E+0 0 \n";
			asciiFile2 << "ellipse(" << map.l(x_b[k], y_b[k]) << "," << map.b(x_b[k], y_b[k]) - (shiftnorth_b?binsize:0.0) << "," << sigma_x[k] << "," << sigma_y[k] << ",0) # color=red width=3\n";
		}
	}
	kold = kold + knew + 10;
	asciiFile1.close();
	asciiFile2.close();
}

ConnectedRegionList* CreateConnectedRegions_ONEBYONE(TH2D* hist, Int_t nconreg, Int_t typeconreg, Int_t nbindeg, Int_t skysegmentation)
{
	ConnectedRegionList* conRegionList = new ConnectedRegionList;
	Int_t M = hist->GetNbinsX(); //dim X
	Int_t N = hist->GetNbinsY();
	conRegionList->M = M;
	conRegionList->N = N;
	Int_t max = (Int_t) hist->GetMaximum() + 1; //si arrotonda
	cout << "--------------------------------" << endl;
	cout << "max " << max << endl;

	Double_t* src = new Double_t[M*N];
	Double_t* src_tmp = new Double_t[M*N];
	for(int i=0; i< M; i++)
		for(int j=0; j< N; j++)
			src[i*M+j] = 0;

	//1) calcola il numero di fotoni di background per ogni pixel
	
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			src_tmp[i * M + j] = hist->GetBinContent(i+1, j+1);

	conRegionList->SetOriginalImage(src_tmp, hist);

	//1.2) e ora calcola il valor medio
	Double_t mean_bkg = 0, points_bkg = 0;
	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src_tmp[i * M + j] != 0)
			{
				points_bkg++;
				mean_bkg += src_tmp[i * M + j];
			}
	mean_bkg /= points_bkg;
	cout << "mean value of the background " << mean_bkg << endl;


	SetImage(M, N, src);
	Int_t last_nconreg = 0;
	int ncon = 0; //number of connected regions


	Int_t dimall = 0;
	for(int k=max; k>=0; k--)
	{
		//salvo l'immagine in src
		for(int i=1; i<= M; i++)
			for(int j=1; j<=N; j++)
			{
				if( ((Int_t)hist->GetBinContent(i,j)) == k)
				{
					src[(i-1)*M + j-1] = hist->GetBinContent(i,j);
					dimall++;
				}
			}
		//per ogni k, provo a calcolare il numero di connessioni connesse

		ncon = LabelingDuePassi(); //in src un histo con le regioni connesse numerate
// 		cout << k << " " << dimall << " " << ncon << endl;

// 		cout << " [" << k << ", " << ncon << "] " ;
		for(int nreg = 1; nreg<=ncon; nreg++)
		{
			Double_t ncon_sig = 0, ncon_pixel = 0;
			for(int i=0; i< M; i++)
				for(int j=0; j<N; j++)
					if(src[i*M + j] == nreg)
					{
						src_tmp[i * M + j] = hist->GetBinContent(i+1, j+1);
						ncon_pixel++;
						ncon_sig += src_tmp[i * M + j];
					}
					else
						src_tmp[i * M + j] = 0;
	
			//si calcola il numero di pixel e l'integrale
			Double_t noise = ncon_pixel * mean_bkg;
			Double_t SNR = ncon_sig / TMath::Sqrt(ncon_sig + 2 * noise);
// 			cout << "************ connected region " << nreg << endl;
// 			cout << "SNR = " << SNR << endl;
// 			cout << "pixel = " << ncon_pixel << endl;
// 			cout << "integral of noise = " << noise << endl;
// 			cout << "integral of signal = " << ncon_sig << endl;
			
		}

		if(ncon >= nconreg)
			break; 

	}

	
	TH2D *bori3 = (TH2D*)hist->Clone("end step 1");
	TString histname = "ESSS1";
	histname += skysegmentation;
	bori3->SetName(histname);

	for(int i=0; i< M; i++)
		for(int j=0; j<N; j++)
			if(src[i * M + j] > 0)
				bori3->SetBinContent(i+1, j+1, src[i * M + j]); //hist->GetBinContent(i+1, j+1)
			else
				bori3->SetBinContent(i+1, j+1, 0);
	bori3->SetTitle("con reg found - number");
	
/*	new TCanvas("con reg found number", "con reg found number");
	bori3->Draw("LEGO2");
	new TCanvas("con reg found number2", "con reg found number2");*/
	
	

	//calcolo della dimensione di ogni regione connessa
	conRegionList->ncon = ncon;
	conRegionList->SetConnectedRegionsImage(src);

	conRegionList->CalculateBaricenterAndError(nbindeg);

	//si fondono le regioni connesse secondo la 8-connection dove una delle due regioni è grande 1
	conRegionList->MergeRegions8Connection();

 	

	conRegionList->CalculateBaricenterAndError(nbindeg);

	
	
	conRegionList->Print();

	nconnectedregion_found = conRegionList->ncon;
	return conRegionList;
}

/*skysegmentation
	0 -> all the sky
	1 -> b > 10
	2 ->  -10 <= b <= 10
	3 -> b < -10
*/


extern void PlotCts2D_FINE(
	bool        shiftnorth_b,
	int         skysegmentation,
	const char* inputfile,
	float       binsize,
	double      sourceradiusremove,
	int         smoothing,
	int         NCONREGSEARCH,
	const char* outFile,
	double      nearradious,
	const char* expfile,
	double      minExp,
	bool        conregsearchType,
	int         removeStep,
	int         nstep,
	int         bkgint,
	const char* filesubtract,
	double      maxMapValue)

{
	ResetArrays();
	/// gStyle->SetPalette(1);	zzz ‘gStyle’ was not declared in this scope
	Int_t nbindeg = 1/binsize;
// 	cout << "nbindeg " << nbindeg << endl;

	AlikeMap map(inputfile);
	double lcenter = map.GetMapCenterL();
	double bcenter = map.GetMapCenterB();

	cerr << "inputfile=" << inputfile << ", lcenter=" << lcenter << ", bcenter=" << bcenter << endl;
	TString histname = "HIST";
	histname += skysegmentation;
	TH2D hist = MakeTH2D(map, histname);
	hist.SetName(histname);
	TH2D hist2(hist);
	histname = "HIST2";
	histname += skysegmentation;
	hist2.SetName(histname);

	if(SetSkyRegion(&hist, skysegmentation, map, sourceradiusremove, lcenter, bcenter) == false)
		return;
	if(SetSkyRegion(&hist2, skysegmentation, map, sourceradiusremove, lcenter, bcenter) == false)
		return;

	int M = hist.GetNbinsX(); //dim X
	int N = hist.GetNbinsY();
	//l'histo 2D come numerazione bin va da (1,1) a (M,N)
// 	cout << "M*N = " << M << " * " << N << endl;

	if (filesubtract[0]) {
		AlikeMap mapSub(filesubtract);

		TH2D histSub = MakeTH2D(mapSub); /// zzz was: TH2D histSub(tFilesubtract);

		histname = "HIST3";
		histname += skysegmentation;
		histSub.SetName(histname);
		if (SetSkyRegion(&histSub, skysegmentation, mapSub, sourceradiusremove, lcenter, bcenter) == false)
			return;
		for(int i=0; i< M; i++)
			for(int j=0; j<N; j++) {
				Float_t sub = TMath::Abs(histSub.GetBinContent(i+1, j+1) - histSub.GetBinContent(i+1, j+1));
				sub = (sub > 0 ? sub : 0);
				hist.SetBinContent(i+1, j+1, sub );
			}
	}

	if(maxMapValue != -1)
		for (int i=0; i< M; i++)
			for (int j=0; j<N; j++)
				if (hist.GetBinContent(i+1, j+1) > maxMapValue) {
// 					cout << "A " << hist->GetBinContent(i+1, j+1) << endl;
					hist.SetBinContent(i+1, j+1, 0);
					hist2.SetBinContent(i+1, j+1, 0);
					}


	if (expfile[0]) {
		AlikeMap imgexp(expfile);
		TH2D histexp = MakeTH2D(imgexp);
		for(int i=1; i<= M; i++)
			for(int j=1; j<=N; j++)
				if(histexp.GetBinContent(i, j) < minExp)
					hist.SetBinContent(i, j, 0);
		}

	hist2.SetTitle("Original");
	#ifdef DRAW
	new TCanvas("Original", "Original");
	hist2.Draw("COL"); //NON CANCELLARE!!!!
	#endif

	if (smoothing)
		Smooth(&hist, M, N, smoothing);

	TH2D bori(hist);
	histname = "HISTBORI";
	histname += skysegmentation;
	bori.SetName(histname);

	bori.SetTitle("Smoothed");
	#ifdef DRAW
	new TCanvas("Smoothed", "Smoothed");
	bori.Draw("COL");
	#endif

	for (int i=0; i<removeStep; i++) {
		cout << "RemoveNoise" << endl;
		RemoveNoise(nstep, &hist, M, N, bkgint);
		}
// 	Int_t dimY;
// 	Int_t dimX;
	hist.Scale(100/hist.GetMaximum()); //scale da 0 a 100
	ConnectedRegionList* regionFinal = new ConnectedRegionList;
// 	TH2D* hh2 = CreateConnectedRegions_ONEBYONE(hist, NCONREGSEARCH, conregsearchType, nbindeg);
// 	AstroImage2 * disp = new AstroImage2(img);
	for (int j=1; j<=NCONREGSEARCH; j++) {
		cout << "STEP " << j << endl;
		ConnectedRegionList* conRegionList = CreateConnectedRegions_ONEBYONE(&hist, j, conregsearchType, nbindeg, skysegmentation);
		for(int k=0; k<MAXLABELCONNECTED; k++) {
			if (conRegionList->list[k].dim != 0) {	
				cout << map.l(conRegionList->list[k].center_x, conRegionList->list[k].center_y) << " " << map.b(conRegionList->list[k].center_x, conRegionList->list[k].center_y) << endl;
			}
		}
		regionFinal->MergeList(conRegionList);
// 		cout << "current version of final list " << endl;
// 		for(int k=0; k<MAXLABELCONNECTED; k++) {
// 			if(regionFinal->list[k].dim != 0) {
// 				cout << disp->Get_l(regionFinal->list[k].center_x, regionFinal->list[k].center_y) << " " << disp->Get_b(regionFinal->list[k].center_x, regionFinal->list[k].center_y) << endl;
// 			}
// 		}
// 		cout << "############### final list " << endl;
// 		regionFinal->Print();
	}
	
	if (nearradious != 0) {
		cout << "Remove near points" << endl;
		regionFinal->RemoveNearPoints(nearradious / binsize);
	}

	cout << "------------------------------ " << endl;
	
	// 	disp->StartPaletteEditor();
	std::ofstream asciiFile1;
	std::ofstream asciiFile2;
	
	TString listFile(outFile);
	//listFile += ".list";
	TString regFile(outFile);
	regFile += ".reg";
	if (skysegmentation == 0 || skysegmentation == 1) {
		asciiFile1.open(listFile);
		asciiFile2.open(regFile);
	} else {
		asciiFile1.open(listFile, std::ios::app);
		asciiFile2.open(regFile, std::ios::app);
	}
	asciiFile2 << "global color=green font=\"helvetica 9 normal\" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\ngalactic\n";

	char outl[50];
	char outb[50];
	cout << "NEW" << endl;

// 	Double_t lcenter = disp->Get_l(bincenter_x, bincenter_y);
// 	Double_t bcenter = disp->Get_b(bincenter_x, bincenter_y) - (shiftnorth_b?binsize:0.0);
	
	cout << "CENTER OF THE MAP " << lcenter << " " << bcenter << endl;	

	static int kold=0;
	int knew = 0;
	for(int k=0; k<MAXLABELCONNECTED; k++) {
		if(regionFinal->list[k].dim != 0) {
			knew = k;
			double lcoo = map.l(regionFinal->list[k].center_x-1, regionFinal->list[k].center_y-1);
			double bcoo = map.b(regionFinal->list[k].center_x-1, regionFinal->list[k].center_y-1) - (shiftnorth_b?binsize:0.0) ;
			sprintf(outl, "%3.2f", lcoo);
			sprintf(outb, "%3.2f", bcoo);
			Double_t dista = AlikeSphdistDeg(lcenter, bcenter, lcoo, bcoo);
			cout << outl << " " << outb << " " << dista << " " << sourceradiusremove<< endl;
			if(sourceradiusremove == 0 || dista <= sourceradiusremove) {
				
				asciiFile1 << "10e-07 " << outl << " " << outb << " 2.1 1 2 S" <<  skysegmentation << k+kold << "\n";
				asciiFile2 << "ellipse(" << outl << "," << outb << "," << regionFinal->list[k].sigma_x << "," << regionFinal->list[k].sigma_y << ",0) # color=red width=1 text={S" << skysegmentation << k+kold << "}\n";
			}
		}
		
	}
	kold = kold + knew + 10;
	asciiFile1.close();
	asciiFile2.close();
}

