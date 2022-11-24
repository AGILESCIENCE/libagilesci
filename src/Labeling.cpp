////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       Author: Andrea Bulgarelli (INAF/IASF Bologna)
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

#include "Labeling.h"

static int     dimX;
static int     dimY;
static double* space_density;

void SetImage(int dim_x, int dim_y, double* space_density_1)
{
	dimX = dim_x;
	dimY = dim_y;
	space_density = space_density_1;
}

static void SetSpaceDensity2D(int i, int j, double value)
{
	if(i < 0 || i >= dimX)
		return;
	if(j < 0 || j >= dimY)
		return;
	space_density[dimY*i+j] = value;
}

static double GetSpaceDensity2D(int i, int j)
{
	if(i < 0 || i >= dimX)
		return 0;
	if(j < 0 || j >= dimY)
		return 0;
	return space_density[dimY*i+j];
}


static int SecondStep(unsigned int last_i, unsigned int NewLabel, int* C )
{
	unsigned int i, k;
	int j, NumGroup;
	int N = dimX;
	int M = dimY;
	//Compattazione degli identificatori delle classi di equivalenza
	//Passaggio non necessario, ma consente una visualizzazione migliore
	NumGroup = 1;
	for(k=1; k<=NewLabel; k++)
	{
		int group = C[k];

		if(group >= NumGroup)
		{
			for(i=k; i<=NewLabel;i++)
				if(C[i] == group) C[i] = NumGroup;
			NumGroup++;
		}
	}

	//secondo passo
	for(i = 0; i < last_i; i++)
		for(j = 0; j < M ; j++)
		{
			//ad ogni pixel con una etichetta gli si assegna
			//una nuova etichetta che coincide con il
			//numero del gruppo (classe di equivalenza) alla
			//quale l'etichetta attuale appartiene
			//			if ((frame[index] != BACKGROUND) && (frame[index] != OBJECT))
			SetSpaceDensity2D(i, j, C[(int)GetSpaceDensity2D(i, j)]);
		}


	for(k=0; k<MAXLABELCONNECTED; k++) C[k] = k;

	return NumGroup;

}

int LabelingDuePassi()
{

	//Versione con soluzione del problema dell'
	//esaurimento delle label

	//IPOTESI DI UTILIZZO:
	// non devono essere presenti oggetti sul bordo;
	// l'immagine deve essere binarizzata.
	//NOTAZIONE
	//			p
	//		q 	x
	//
	//(4-connettivita')

	int lx = 0;		//label del pixel x
	unsigned int NewLabel = 0;	//tiene traccia delle label non ancora utilizzate
	int C[MAXLABELCONNECTED];		//Array di equivalenza
	int lp, lq;		//label dei pixel p e q
	unsigned int last_i;

	int i, j;
	unsigned int k;
	int NumObject = 0;	//utilizzato per contare il numero di oggetti individuati
	// 	UChar_t BG = 0;
	long int count=0;
	long int cc = 0;
	long int bb = 0;

	int N = dimX;
	int M = dimY;

	int F[MAXLABELCONNECTED];

	//init di R
	for(k=0; k<MAXLABELCONNECTED; k++) C[k] = k, F[k] = 1;

	//prima fase
	for(i=0; i<N; i++)
		for(j=0; j<M; j++)
		{
			if (GetSpaceDensity2D(i, j) != BG)	//e' un punto da analizzare
			{

				bb++;
				lp = GetSpaceDensity2D(i-1, j);
				lq = GetSpaceDensity2D(i, j-1);
				cc++;
				if(C[lp] == C[lq])
				{
					cc++;
					if(lp != BG) lx = lp;
					else
					{	//creazione di una nuova label
						NewLabel++;
						NumObject++;

						if(NewLabel >= 255)
						{
							NewLabel --;
							last_i = i;
							NewLabel = SecondStep(last_i, NewLabel, C);
							for(k=0; k<MAXLABELCONNECTED; k++) F[k] = 1;

							if(NewLabel >= MAXLABELCONNECTED)
							{
								//NotLabeling();
								return -1;
							}
						}

						lx = NewLabel;
					}
				}
				else
				{
					//C[lp] != C[lq]
					cc++;
					if(lp == BG ) lx = lq;
					else
					{
						cc++;
						if(lq == BG)
							lx = lp;
						else	//FONDI
						{
							int Clp = C[lp], Clq = C[lq];

							//unione delle due classi di equivalenza Clp e Clq
							if(F[Clq] == 1)		//e' presente una sola label
							{	C[lq] = Clp;	//si fonde Clq eliminandolo
								//con una sola sostituzione
								F[Clp] = 0;
								//F[Clq] = 0;
							}
							else if(F[Clp] == 1)
							{
								C[lp] = Clq;

								F[Clp] = 0;
								//F[Clq] = 0;
							}
							else
							{
								for(k=1; k<= NewLabel; k++)
									if(C[k] == Clp)
										C[k] = Clq; //elimina Clp
								F[Clq] = 0; //sono presenti pi label;
								count++;
							}

							lx = lq;	//scelta di una delle due label
							NumObject--;
							//					count++;

						}
					}
				}
				SetSpaceDensity2D(i, j, lx);
			}
		}

	last_i = N ;

	NewLabel = SecondStep(last_i, NewLabel, C);

	return NumObject;

}

int LabelingDuePassi8()
{

	//Versione con soluzione del problema dell'
	//esaurimento delle label

	//NOTAZIONE
	//		p	q	r
	//		s 	x
	//
	//(8-connettivita')

	int lx = 0;		//label del pixel x
	unsigned int NewLabel = 0;	//tiene traccia delle label non ancora utilizzate
	int C[MAXLABELCONNECTED];		//Array di equivalenza
	int lp, lq, lr, ls;		//label dei pixel p e q
	unsigned int last_i;

	int i, j;
	unsigned int k;
	int NumObject = 0;	//utilizzato per contare il numero di oggetti individuati
	// 	UChar_t BG = 0;
	long int count=0;

	int N = dimX;
	int M = dimY;

	int F[MAXLABELCONNECTED];

	//init di R
	for(k=0; k<MAXLABELCONNECTED; k++) C[k] = k, F[k] = 1;

	//prima fase
	for(i=0; i<N; i++)
		for(j=0; j<M; j++)
		{
			if (GetSpaceDensity2D(i, j) != BG)	//e' un punto da analizzare
			{

				lp = GetSpaceDensity2D(i-1, j-1);
				lq = GetSpaceDensity2D(i-1, j);
				lr = GetSpaceDensity2D(i-1, j+1);
				ls = GetSpaceDensity2D(i, j-1);

				if(C[lp] == C[lq] && C[lq] == C[lr] && C[lr] == C[ls])
				{
					if(lp != BG)
						lx = ls;
					else
					{
						//tutte le classi di equivalenza sono uguali, ma sono BG
						//creazione di una nuova label
						NewLabel++;
						NumObject++;

						if(NewLabel >= 255)
						{
							NewLabel --;
							last_i = i;
							NewLabel = SecondStep(last_i, NewLabel, C);
							for(k=0; k<MAXLABELCONNECTED; k++) F[k] = 1;

							if(NewLabel >= MAXLABELCONNECTED)
							{
								//NotLabeling();
								return -1;
							}
						}

						lx = NewLabel;
					}
				}
				else
				{
					//almeno una delle classi di equivalenza è diversa.
					//NB: se si è arrivati fino a qui le classi di equivalenza
					//di p, q ed r devono essere la stessa, a meno che p o q o r
					//non appartengano al BG, oppure non si presenti il seguente caso (A)
					//	p	BG	r
					//	BG	x
					// Si sfrutta questa caratteristica
					//A) Si verifica se siamo nel caso A
					if(C[lp] != C[lr] && lp != BG && lr!= BG)
					{
						//esiste x, allora fondi
						int Clp = C[lp], Clr = C[lr];

						//unione delle due classi di equivalenza Clp e Clq
						if(F[Clr] == 1)		//e' presente una sola label
						{	C[lr] = Clp;	//si fonde Clq eliminandolo
							//con una sola sostituzione
							F[Clp] = 0;
							//F[Clq] = 0;
						}
						else if(F[Clp] == 1)
						{
							C[lp] = Clr;

							F[Clp] = 0;
							//F[Clq] = 0;
						}
						else
						{
							for(k=1; k<= NewLabel; k++)
								if(C[k] == Clp)
									C[k] = Clr; //elimina Clp
							F[Clr] = 0; //sono presenti pi label;
							count++;
						}
						NumObject--;
						lx = lr;
						//					count++;
					}
					else
					{

						//B) Casp B. Si cerca la prima != da BG tra lp, lq e lr
						int lk = BG;
						if(lp != BG)
							lk = lp;
						if(lq != BG)
							lk = lq;
						if(lr != BG)
							lk = lr;
						//si confronta lk con ls
						if(lk == BG )
							lx = ls;
						else
							if(ls == BG)
								lx = lk;
							else	//FONDI
							{
								int Clk = C[lk], Cls = C[ls];

								//unione delle due classi di equivalenza Clp e Clq
								if(F[Cls] == 1)		//e' presente una sola label
								{	C[ls] = Clk;	//si fonde Clq eliminandolo
									//con una sola sostituzione
									F[Clk] = 0;
									//F[Clq] = 0;
								}
								else if(F[Clk] == 1)
								{
									C[lk] = Cls;

									F[Clk] = 0;
									//F[Clq] = 0;
								}
								else
								{
									for(k=1; k<= NewLabel; k++)
										if(C[k] == Clk)
											C[k] = Cls; //elimina Clk
									F[Cls] = 0; //sono presenti pi label;
									count++;
								}

								lx = ls;	//scelta di una delle due label
								NumObject--;
								//					count++;

							}
					}
				}
				SetSpaceDensity2D(i, j, lx);
			}
		}

	last_i = N ;

	NewLabel = SecondStep(last_i, NewLabel, C);

	return NumObject;

}
