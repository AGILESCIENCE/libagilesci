////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       Authors: Leonardo Baroncelli, Giancarlo Zollino, 2017
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


#include "BinEvaluator.h"

BinEvaluator::BinEvaluator(const char *_fitsFilePath, DoubleMatrixCustomMap * _image, double _l, double _b, double _radius) {
	image = _image;
	fitsFilePath=_fitsFilePath;
	l=_l;
	b=_b;
	radius=_radius;
	binSum=0;
	agileMapUtils = new AgileMap(_fitsFilePath);
	tmin = agileMapUtils->GetTstart();
	tmax = agileMapUtils->GetTstop();
	x=0;
	y=0;
	agileMapUtils->GetRowCol(l,b,&x,&y);
	rows = agileMapUtils->Rows();
	cols = agileMapUtils->Cols();
}

BinEvaluator::BinEvaluator(const char * _fitsFilePath, double _l, double _b, double _radius) {
	image = MapConverter::fitsMapToDoubleMatrix(_fitsFilePath);
	fitsFilePath=_fitsFilePath;
	l=_l;
	b=_b;
	radius=_radius;
	binSum=0;
	agileMapUtils = new AgileMap(_fitsFilePath);
	tmin = agileMapUtils->GetTstart();
	tmax = agileMapUtils->GetTstop();
	x=0;
	y=0;
	agileMapUtils->GetRowCol(l,b,&x,&y);
	rows = agileMapUtils->Rows();
	cols = agileMapUtils->Cols();

}

/* No need to delete image, because it is deleted by ExpRatioEvaluator destructor.
BinEvaluator::~BinEvaluator(){
	delete agileMapUtils;
}*/

void BinEvaluator::sumBin()
{
	double greyLevel;
	binSum = 0;
	if(checkIfRadiusExceedsMapBorders()){
		for(int i = 0; i < rows; i++){
			for(int j=0; j < cols; j++){
				greyLevel = image->image[i][j];
				if( greyLevel > 0  &&  agileMapUtils->SrcDist(i,j,l,b) < radius){
					binSum += greyLevel;
				}
			}
		}
	}
	else{
		fprintf(stderr,"[BinEvaluator Error]: %s -> %f the radius exceeds the border of the .exp map\n",fitsFilePath,radius);
		exit(EXIT_FAILURE);
	}
}



bool BinEvaluator::checkIfRadiusExceedsMapBorders() {

	double distSx;
	double distDx;
	double distUp;
	double distDown;

	distSx   = sqrt(pow(double(0-x),2));
	distDx   = sqrt(pow(double(cols-1-x),2));
	distUp   = sqrt(pow(double(0-y),2));
	distDown = sqrt(pow(double(rows-1-y),2));

	if( distSx < radius || distDx < radius || distUp < radius || distDown < radius )
		return false;
	else
		return true;
 }
