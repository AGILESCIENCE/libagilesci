
/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
 *
*/

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
