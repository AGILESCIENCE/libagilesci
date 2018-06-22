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

#include <string>
#include <sstream>
#include <iomanip>
#include "AgileMap.h"
#include "FitsUtils.h"
#include "AlikeData5.h"
#include "MapConverter.h"

using std::string;
using std::ofstream;
using std::ios;
using std::setprecision;

class ExpRatioEvaluator
{
	public:

	////////////////////////////////////////////////////////////////////////
	//
	// 	PUBLIC CONSTRUCTORS
	//
	//	If onNormalizedMap == true -> creates normalizedImage and writes it on file .exp.gz
	//	If createExpRatioMap == true -> creates createExpRatioMap and writes it on file .exp.gz

		// Creates image extracting pixel values from expPath.
		ExpRatioEvaluator(const char * expPath,bool isExpMapNormalized, bool createExpNormalizedMap, bool createExpRatioMap, double minThreshold, double maxThreshold, double squareSize);

		// Creates image extracting pixel values from AgileMap.
		ExpRatioEvaluator(AgileMap agileMap, bool isExpMapNormalized, bool createExpNormalizedMap, bool createExpRatioMap, double minThreshold, double maxThreshold, double squareSize);

		//
		~ExpRatioEvaluator();




	//////////////////////////////////////////////////////////////////////////
	//
	// 	PUBLIC METHODS
	//
	//
	//  Computes and returns the output array. If the rectangle is not entirely inside the image, it returns -1 -1 -1 -1.
		double computeExpRatioValues(double l, double b);
		double computeExpRatioValues(int x, int y, string type);





	// Getters
		DoubleMatrixCustomMap * getImage();
		DoubleMatrixCustomMap * getExpRatioMap();
		DoubleMatrixCustomMap * getNormalizedMap();

	// Get filename
		const char * getExpPath();

	// Get ExpRatio Parameters
		double getMinThreshold();
		double getMaxThreshold();
		double getSquareSize();





	private:

		// Utility class
		AgileMap agileMap;

		const char* expPath;

		// Degrees pixel size
		double cdelt2;

		// If pixel value is not inside [minThreshold, maxThreshold], increments nBad.
		double minThreshold;
		double maxThreshold;

		// The size of the rectangle (x-size , x+size, y-size, y+size)
		double squareSize;

		// If false, the exp map given in input must be normalized
		bool isExpMapNormalized;

		// If true, creates the expRatio map
		bool createExpRatioMap;
		bool createExpNormalizedMap;


		// The images
		DoubleMatrixCustomMap * image;
		DoubleMatrixCustomMap * normalizedImage;
		DoubleMatrixCustomMap * expRatioImage;
		
		DoubleMatrixCustomMap * normalizationFactorMatrix;

		// The output array  [ exp-ratio, nBad, nTot, greyLevelMean ]
		double output[4];




		void createAndWriteImages();


		// Check if the  rectangle is completely inside the image
		bool isRectangleInside(int x, int y);

		// Computes and returns the output array. If the rectangle is not entirely inside the image, it returns -1 -1 -1 -1.
		double computeExpRatio(int x, int y);


		// Used function to compute normalizationFactorMatrix
		double Alikesinaa(double input);



	////////////////////////////////////////////////////////////////////////
	//
	// Initialize a distance matrix.
	// Computes a distance matrix from double ** image using SrcDist( i, j, center_l, center_b )
	// to store the distances from every pixel of the first quadrant of the matrix (i,j)
	// to the center (center_l, center_b).
	//
	// Initialize and compute the normalizationFactorMatrix, iterating ONLY ON THE FIRST QUADRANT
	// OF THE MATRIX, WRITING 4 PIXEL (the symmetric ones) AT ONE TIME, with this formula:
	//
	// if the distance is not 0
	// 0.0003046174197867085688996857673060958405 * xbin * ybin * Alikesinaa(0.0174532925199432954743716805978692718782 * distanceMatrix[i][j]);
	// else
	// 0.0003046174197867085688996857673060958405 * xbin * ybin

		DoubleMatrixCustomMap * computeSpatialNormalizationFactorMatrix();


	// Compute the normalizedImage using this formula:
	// normalizedImage[i][j] = normalizedImage[i][j] / timeNormalizationFactor * normalizationFactorMatrix[i][j]
		DoubleMatrixCustomMap * createNormalizedImage();


	// Computes expRatioImage in which the value of each pixel is the exp-ratio value centered in the same pixel.
		DoubleMatrixCustomMap * createExpRatioPixelMap();




	// Writes on disk a fits file.
		int writeMatrixDataInAgileMapFile(const char * appendToFileName, DoubleMatrixCustomMap * matrixData);
		string computeNewFileName(const char * appendToFilename);
};
