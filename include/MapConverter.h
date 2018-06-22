/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino,
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 *
 * https://github.com/Leofaber/MapConverter
*/

#ifndef MAPCONVERTER_H
#define MAPCONVERTER_H

#include "FitsUtils.h"
#include <iostream>
#include <string>

using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ostringstream;

class IntMatrixCustomMap {

	public:
		int ** image;
		int rows;
		int cols;
		string imagePath;

		IntMatrixCustomMap(){
			rows = 0;
			cols = 0;
			imagePath = "";
		}


		~IntMatrixCustomMap(){
			for (  int i = 0; i < rows; i++){
				delete [] image[i];
			}
			delete [] image;
			image = 0;
		}


		void initializeImage(string imgPath, int r, int c){
			imagePath = imgPath;
			rows = r;
			cols = c;
			image = new int*[rows];
			for (int y = 0; y < rows; ++y){
				image[y] = new int[cols];
			}
		}

};


class DoubleMatrixCustomMap {

	public:
		double ** image;
		int rows;
		int cols;
		string imagePath;

		DoubleMatrixCustomMap(){
			rows = 0;
			cols = 0;
			imagePath = "";
		}


		~DoubleMatrixCustomMap(){
			for (  int i = 0; i < rows; i++){
				delete [] image[i];
			}
			delete [] image;
			image = 0;
		}


		void initializeImage(string imgPath, int r, int c){
			imagePath = imgPath;
			rows = r;
			cols = c;
			image = new double*[rows];
			for (int y = 0; y < rows; ++y){
				image[y] = new double[cols];
			}
		}

};


class MapConverter
{
	public:


		// Convert an image in a **double matrix
		static DoubleMatrixCustomMap * fitsMapToDoubleMatrix(const char * fitsImagePath);

		// Convert an image in a **double matrix
		static IntMatrixCustomMap * fitsMapToIntMatrix(const char * fitsImagePath);



	private:
		// Constructor
		MapConverter();
};

#endif
