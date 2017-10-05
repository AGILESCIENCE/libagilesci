/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino,
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 * 
 * https://github.com/Leofaber/ExpRatioEvaluator
*/

#include "ExpRatioEvaluator.h"


ExpRatioEvaluator::ExpRatioEvaluator(const char * _expPath) 
{
	expPath=_expPath;
	
	agileMap=new AgileMap(expPath);
	
	tStart=agileMap->GetTstart();
	tStop=agileMap->GetTstop();
	timeFactor=tStop-tStart;
 
	double cdelt2=agileMap->GetYbin();
	size = 10/cdelt2;
	
	
	spatialFactor = 0.0003046174197867085688996857673060958405*cdelt2*cdelt2;
	normalizationFactor = spatialFactor*timeFactor;
	
	// Initialization of image and normalizedImage
	if(! convertFitsDataToMatrix() )
	{
		fprintf( stderr, "[ExpRatioEvaluator] ERROR!! convertFitsDataToMatrix(): error reading fits file\n");
		exit (EXIT_FAILURE);
	}


	

	// Writing of normalizedImage
	if( writeMatrixDataInAgileMapFile( normalizedImage, agileMap,  "norm.exp") )	
		cout << "\nCreated AgileMap file: exp_norm.exp" << endl;
	else
		cout << "\nCAN'T create AgileMap file: exp_norm.exp" << endl;
	
} 


bool ExpRatioEvaluator::convertFitsDataToMatrix()
{
	


	//CFITSIO
	fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	int bitpix, naxis, ii, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	double *pixels;
	char format[20], hdformat[20];
			

	if (!fits_open_file(&fptr, expPath, READONLY, &status))
	{									// 16   , 2     , {166,166}
		if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status))
		{
			if (naxis > 2 || naxis == 0)
			{
				printf("Error: only 1D or 2D images are supported\n");
				return false;
			}			
			else
			{	 
				rows = (int)naxes[0]; 
				cols = (int)naxes[1];
				image = new double*[rows];
 				normalizedImage = new double*[rows];

				for (int i = 0; i < rows; ++i){
					image[i] = new double[cols];
					normalizedImage[i] = new double[cols];
				}

				/* get memory for 1 row */
				pixels = (double *)malloc(naxes[0] * sizeof(double));

				if (pixels == NULL)
				{
					printf("Memory allocation error - Press any key to exit\n");
					return false;
				}
				else
				{
					/* loop over all the rows in the image, top to bottom */

					int col_index = 0;
					int row_index = 0;
					for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
					{
						if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status))  /* read row of pixels */
							break;  /* jump out of loop on error */

						for (ii = 0; ii < naxes[0]; ii++)
						{
							
							image[row_index][col_index] = (double)pixels[ii];
							normalizedImage[row_index][col_index] = (double)pixels[ii]/normalizationFactor;
							//cout << image[row_index][col_index] << " ";
							col_index++;
						}
						col_index = 0;
						row_index++;
						//cout << "\n";
					}

					free(pixels);
				}
			}

		}

		fits_close_file(fptr, &status);

	}
	if (status>0)
	{
		printf("Can't open fits file - Press any key to exit\n");
		return false;	
	}

	return true;

}


bool ExpRatioEvaluator::isRectangleInside(int x, int y) 
{
	double distSx;
	double distDx;
	double distUp;
	double distDown;

	distSx =  pow(pow(double(0-x),2),0.5);
	distDx =  pow(pow(double(cols-1-x),2),0.5);
	distUp =  pow(pow(double(0-y),2),0.5);
	distDown =  pow(pow(double(rows-1-y),2),0.5);
	if(distSx < size || distDx < size || distUp < size || distDown < size)
		return false;
	else
		return true;

	
}




double* ExpRatioEvaluator::computeExpRatioValues(double l, double b, bool onNormalizeMap, double minThreshold, double maxThreshold) 
{ 
	int x;
	int y;

	agileMap->GetRowCol(l,b,&x,&y);
	if(x < 0 || x > cols-1 || y < 0 || y > rows-1)
	{
		fprintf( stderr, "[ExpRatioEvaluator] ERROR!! computeExpRatioValues(double l, double b, bool onNormalizeMap, double minThreshold, double maxThreshold):  Map .cts and Map .exp are not centered in the same galactic coordinate because input l and input b go out of the .exp map\n");
		output[0] =  -2;
		output[1] =  -2;
		output[2] =  -2;
		output[3] =  -2;
		return output;
	}
	return computeExpRatio(x, y, onNormalizeMap, minThreshold, maxThreshold);
}


double* ExpRatioEvaluator::computeExpRatioValues(int x, int y, string type, bool onNormalizeMap, double minThreshold, double maxThreshold)
{	
	return computeExpRatio(x, y, onNormalizeMap, minThreshold, maxThreshold);
}


double * ExpRatioEvaluator::computeExpRatio(int x, int y, bool onNormalizeMap, double minThreshold, double maxThreshold){
	
	/*cout << "MinThreshold: "<< minThreshold << endl;
	cout << "MaxThreshold: " << maxThreshold << endl;
	if(onNormalizeMap)
		cout << "ExpRatioEvaluator will compute exp-ratio value of the normalized map" << endl;
	else	
		cout << "ExpRatioEvaluator will compute exp-ratio value of the NO-normalized map" << endl;	*/

 	int xmin, xmax, ymin,ymax;
	int npixel = 0;
	int nBad = 0;
	double totCount=0;
	double tmp = 0;
	double expRatio = 0;
	double greyLevelSum = 0;
	double mean = 0;
		
	xmin = x - size;
	xmax = x + size;
	ymin = y - size;
	ymax = y + size;


	if(isRectangleInside(x,y)) 
	{	
		for(int i = xmin; i <= xmax; i++) 
		{
			for(int j= ymin; j <= ymax; j++) 
			{	
				totCount+=1;
				if(onNormalizeMap) 
					tmp=(double)normalizedImage[i][j];
				else
					tmp=(double)image[i][j];

				greyLevelSum+=tmp;

				if(tmp < minThreshold || tmp >maxThreshold)
					nBad+=1;
			}
		}
	
		output[0] = (1-(nBad/totCount))*100;
		output[1] = nBad;
		output[2] = totCount;
		output[3] = greyLevelSum/totCount;
		return output;

	}else 
	{
		//fprintf( stderr, "Rectangle is not completely inside!\n");
		output[0] =  -1;
		output[1] =  -1;
		output[2] =  -1;
		output[3] =  -1;
		return output;

	}
}


 
int ExpRatioEvaluator::getRows(){
	return rows;
}
int ExpRatioEvaluator::getCols(){
	return cols;
}
																													// _norm.exp
bool ExpRatioEvaluator::writeMatrixDataInAgileMapFile(double ** matrixData, AgileMap * agileMapForCopy, const char * appendToFilename){
	// copio il file AgileMap
	AgileMap* newMap = new AgileMap(*agileMapForCopy);
	
	const char * imageName = agileMapForCopy->GetFileName(); // e.g.    MAP1000s.exp
	string imageName_string(imageName);
	
	string newFileName = "";
	
    size_t foundPatternExp = imageName_string.find(".exp");
	size_t foundPatternCts = imageName_string.find(".cts");

    if(foundPatternExp != string::npos)
		newFileName = imageName_string.substr(0,foundPatternExp);
	else if(foundPatternCts != string::npos)
		newFileName = imageName_string.substr(0,foundPatternCts);
    else
		newFileName = imageName_string;

	string appendToFilenameString(appendToFilename);
	newFileName +="_"+appendToFilenameString;

		
	const char * newFileNameC = newFileName.c_str();
	remove( newFileNameC );

	// lo scrivo su file cambiandogli nome
	int statusWrite = newMap->Write(newFileNameC);

	// copio i dati nel file aprendolo con cfitsio
	bool statusCopy = copyDataToFitsFile(newFileNameC, matrixData);

	if(statusWrite == 0 && statusCopy)
		return true;
	else
		return false;

}



bool ExpRatioEvaluator::copyDataToFitsFile(const char * pathToAgileMapFile, double ** matrixData){
	
	//CFITSIO
	fitsfile *fptr;    
	int status = 0;    
	int bitpix, naxis, ii, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	double *pixels;
	char format[20], hdformat[20];
			

	if (!fits_open_file(&fptr, pathToAgileMapFile, READWRITE, &status))
	{									// 16   , 2     , {166,166}
		if (!fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status))
		{
			 		
 
			pixels = (double *)malloc(naxes[0] * sizeof(double));

 

			int col_index = 0;
			int row_index = 0;
			for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
			{

//int fits_read_pix (fitsfile *fptr, int datatype, long *fpixel, LONGLONG nelements,DTYPE *nulval, > DTYPE *array, int *anynul, int *status)
//int fits_write_pix(fitsfile *fptr, int datatype, long *fpixel, LONGLONG nelements,DTYPE *array, int *status);

				for (ii = 0; ii < naxes[0]; ii++)
				{	 
				 	pixels[ii] = matrixData[row_index][col_index];
				 	col_index++;
				}
				fits_write_pix(fptr, TDOUBLE, fpixel, naxes[0], pixels ,     &status);
				col_index = 0;
				row_index++;

			}

			free(pixels);		 

		}
		else
		{
			printf("Can't read fits params\n");
			return false;
		}

		fits_close_file(fptr, &status);

	}else
	{
		printf("Can't open fits file\n");
		return false;	
	}

	return true;
	
}

double ** ExpRatioEvaluator::createExpRatioPixelMap(bool computeExpRatioOnNormalizedMap, double minThreshold, double maxThreshold){



		// GENERAZIONE DELLA MAPPA IN CUI OGNI PIXEL E' UN EXP-RATIO
			
		// il seguente codice, data la mappa normalizzata, crea una seconda mappa, in cui ogni pixel è il valore di exp-ratio centrato
		// sul pixel corrispondente della mappa normalizzata. (i bordi saranno -1 perchè il rettangolo fuoriuscirà).
		// La nuova mappa deve essere inserita in una file CTS leggibile con ds9.

		// dovrebbe uscire una roba del genere:

		// -1 -1 -1 -1 -1 -1 ........
		// -1 -1 -1 -1 -1 -1 ........
		// -1 -1... .. .. .. ........
		// -1 -1...  a  b  c ........
		// -1 -1...  q  w  e ........
		// .. .....  x  y  z ........
		// .. ..... .. .. .. ........


		// Algoritmo:
		
		// crea una nuova mappa e inizializzala
		double ** expRatioMap;
		expRatioMap = new double*[rows];
				
		for(int i=0; i<rows; i++) {
			expRatioMap[i] = new double[cols];
		}
					
		
		// calcola i valori e inseriscili nella nuova mappa
		for(int i = 0; i < rows ; i++ ) {
			for(int j = 0; j < cols; j++) {
 
				double *output2 = computeExpRatioValues(i,j,"PIXEL",computeExpRatioOnNormalizedMap,minThreshold,maxThreshold);
				expRatioMap[i][j] = output2[0];
 
			}
			
		}
		
		// stampa (debug)
	/*	for(int i = 0; i<rows ; i++ ) {
			for(int j = 0; j < cols; j++) {
				cout << expRatioMap[i][j] << " ";
			} 
			cout << "\n";
		}*/
		
		return expRatioMap;


}
