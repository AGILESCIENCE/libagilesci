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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>

#include "AgileMap.h"
#include "FitsUtils.h"
#include "AlikeData5.h"

using namespace std;

class ExpRatioEvaluator
{
	public: 
	
		/*
			Creates image extracting pixel values from expPath.
			Creates normalizedImage and expRatioImage, and writes them on file .exp.gz 
		*/
		ExpRatioEvaluator(const char * expPath,bool onNormalizedMap, double minThreshold, double maxThreshold, int squareSize);

		/*
			Constructor used in multi5 project
			Creates image extracting pixel values from expPath.
			Creates normalizedImage and expRatioImage, and writes them on file .exp.gz
		*/
		ExpRatioEvaluator(const char * expPath,bool onNormalizedMap);

		/*
			Computes and returns the output array. 
			If the rectangle is not entirely inside the image, it returns -1 -1 -1 -1.
		*/	
		double computeExpRatioValues(double l, double b);	
		double computeExpRatioValues(int x, int y, string type);

		
		// Getting Map
		double ** getImage();
		double ** getExpRatioMap();
		double ** getNormalizedMap();

		// Get filename
		const char* getExpPath();

		// Get ExpRatio Parameters
		double getMinThreshold();
		double getMaxThreshold();
		int getSquareSize();

	private:

		/*
			**************** ATTRIBUTES ****************
		*/

		// Utility class	
		AgileMap* agileMap;

		const char* expPath;
	
		// Degrees pixel size
		double cdelt2;

		// If pixel value is not inside [minThreshold, maxThreshold], increments nBad.
		double minThreshold;
		double maxThreshold;
	
		// The size of the rectangle (x-size , x+size, y-size, y+size)
		float squareSize;

		// If true, the expratio evaluation is done on the normalizedImage 
		bool onNormalizedMap;

		// Lenght of fits data axes
		int rows;
		int cols;

		// Input for Alikesinaa
		//double input;


		// The images 	
		double ** image;
		double ** normalizedImage;
		double ** expRatioImage;
		double ** normalizationFactorMatrix;

		// The output array  [ exp-ratio, nBad, nTot, greyLevelMean ]	
		double output[4];

		


	
		/*
			**************** METHODS ****************
		*/
	

		// Open the exp image and extract the pixel values and store them in the double** image
		bool convertFitsDataToMatrix();
	
		// Check if the  rectangle is completely inside the image
		bool isRectangleInside(int x, int y);

		/*
			Computes and returns the output array. 
			If the rectangle is not entirely inside the image, it returns -1 -1 -1 -1.
		*/	
		double computeExpRatio(int x, int y);	


		/* 	
			Used function to compute normalizationFactorMatrix
		*/
		double Alikesinaa(double input);



		/*			
			Initialize a distance matrix.

			Computes a distance matrix from double ** image using SrcDist( i, j, center_l, center_b ) to store the distances from every pixel of
			the first quadrant of the matrix (i,j) to the center (center_l, center_b).

			Initialize and compute the normalizationFactorMatrix, iterating ONLY ON THE FIRST QUADRANT OF THE MATRIX, WRITING 4 PIXEL (the symmetric 				ones) AT TIME, with this formula:

			if the distance is not 0 
			0.0003046174197867085688996857673060958405 * xbin * ybin * Alikesinaa(0.0174532925199432954743716805978692718782 * distanceMatrix[i][j]);
			else
			0.0003046174197867085688996857673060958405 * xbin * ybin

		*/	
		double ** computeSpatialNormalizationFactorMatrix();

	
		/*
			Compute the normalizedImage using this formula:
			normalizedImage[i][j] = normalizedImage[i][j] / timeNormalizationFactor * normalizationFactorMatrix[i][j]
		*/
		double ** createNormalizedImage();


		/*
			Computes expRatioImage in which the value of each pixel is the exp-ratio value centered in the same pixel. 
		*/
		double ** createExpRatioPixelMap(double minThreshold, double maxThreshold);

		

	
		/*
			Crea un file AgileMap (apribile con ds9). Se il file originale si chiama MAP.exp, e appendToFileName="norm_exp.exp" 
			il file che verrà scritto su disco avrà nome MAP_norm_exp.exp.gz
		*/
		int writeMatrixDataInAgileMapFile(const char * appendToFileName, double ** matrixData);
	
		string computeNewFileName(const char * appendToFilename);
};
