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


#include <iostream>
#include <stdlib.h> 
#include <stdio.h>
#include "fitsio.h"
#include "AgileMap.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

class ExpRatioEvaluator
{
	public: 

	ExpRatioEvaluator(const char * expPath);


	const char* expPath;
	double normalizationFactor;
	double tStart;
	double tStop;
	double timeFactor;
	double spatialFactor;
	
	
	// If pixel value is not inside [minThreshold, maxThreshold], increments nBad.
	double minThreshold;
	double maxThreshold;
	
	AgileMap* agileMap;

	// The size of the rectangle (x-size , x+size, y-size, y+size)
	float size;

	// Check if the  rectangle is completely inside the image
	bool isRectangleInside(int x, int y);

	// We convert fits data into a matrix of double
	int rows;
	int cols;
	double ** image;
	double ** normalizedImage;

	bool convertFitsDataToMatrix();
	
	// The output array  [ exp-ratio, nBad, nTot, greyLevelMean ]	
	double output[4];



	/*
		Computes and returns the output array. 
		If the rectangle is not entirely inside the image, it returns -1 -1 -1 -1.
	*/	
	double* computeExpRatioValues(double l, double b, bool onNormalizeMap, double minThreshold, double maxThreshold);	
	double* computeExpRatioValues(int x, int y, string type, bool onNormalizeMap, double minThreshold, double maxThreshold);
	
	double* computeExpRatio(int x, int y, bool onNormalizeMap, double minThreshold, double maxThreshold);	

	/*
		Crea un file AgileMap (apribile con ds9). Per far ciò serve un file AgileMap da copiare, una matrice da cui estrarre i dati da inserire nel nuovo file, ed il nome del nuovo file.
	*/
	bool writeMatrixDataInAgileMapFile(double ** matrixData, AgileMap * agileMapForCopy, const char * appendToFilename);

	/*
		Apre con cfitsio il file in posizione pathToAgileMapFile e scrive all'interno i valori estratti da matrixData
	*/
	bool copyDataToFitsFile(const char * pathToAgileMapFile, double ** matrixData);


	/*
		Crea una mappa in cui il valore di ogni pixel è l'exp-ratio calcolato sul pixel stesso. 
	*/
	double ** createExpRatioPixelMap(bool computeExpRatioOnNormalizedMap, double minThreshold, double maxThreshold);


	 
	int getRows();
	int getCols();

	 
	
		
	

};
