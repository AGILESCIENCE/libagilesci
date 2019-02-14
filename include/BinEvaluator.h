////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       2017
//       Authors: Leonardo Baroncelli, Giancarlo Zollino (INAF/OAS Bologna)
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

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <AgileMap.h>
#include <fitsio.h>
#include <stdio.h>
#include <iostream>
#include "MapConverter.h"


class BinEvaluator {
	public:

		BinEvaluator(const char * fitsFilePath, DoubleMatrixCustomMap * image, double l, double b, double radius);
		BinEvaluator(const char * fitsFilePath, double l, double b, double radius);
        //~BinEvaluator();

        double binSum;
        double tmin;
        double tmax;

        void sumBin();

    private:

        const char * fitsFilePath;

        int    x, y;
        double l, b;

        double radius;

        int rows;
        int cols;
        AgileMap * agileMapUtils;

        DoubleMatrixCustomMap * image;

        bool checkIfRadiusExceedsMapBorders();


};

#endif
