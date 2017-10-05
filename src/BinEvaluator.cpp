/*
 * Copyright (c) 2017
 *     Leonardo Baroncelli, Giancarlo Zollino
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/


#include "BinEvaluator.h"

BinEvaluator::BinEvaluator(const char *_fitsFilePath, double _l, double _b, double _radius) {
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
 
}

bool BinEvaluator::convertFitsDataToMatrix() {
	
//CFITSIO
	fitsfile *fptr;   /* FITS file pointer, defined in fitsio.h */
	int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
	int bitpix, naxis, ii, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	double *pixels;
	char format[20], hdformat[20];
			

	if (!fits_open_file(&fptr, fitsFilePath, READONLY, &status))
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
				for (int i = 0; i < rows; ++i){
					image[i] = new double[cols];
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
					{	/* read row of pixels */
						if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status)){
							printf("fits_read_pix error \n");
							return false; 
						}  
							 /* jump out of loop on error */

						for (ii = 0; ii < naxes[0]; ii++)
						{
							image[row_index][col_index] = (double)pixels[ii];
							col_index++;
						}
						col_index = 0;
						row_index++;
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



int BinEvaluator::sumBin() {
	
	
	int status,i,j;
	double greyLevel;

	
	if(isRadiusInside()) { 
		for(int i = 0; i < rows; i++){
			for(int j=0; j < cols; j++){
					greyLevel = image[i][j];
					if(greyLevel>0 && agileMapUtils->SrcDist(i,j,l,b)<=radius){
						binSum+=greyLevel;

				}
			}
		}
	}else{
		return -1;
	}
	return 0;
}

bool BinEvaluator::isRadiusInside() {
	
	double distSx;
	double distDx;
	double distUp;
	double distDown;

	distSx =  sqrt(pow(double(0-x),2));
	distDx =  sqrt(pow(double(cols-1-x),2));
	distUp =  sqrt(pow(double(0-y),2));
	distDown = sqrt(pow(double(rows-1-y),2));
	if(distSx < radius || distDx < radius || distUp < radius || distDown < radius)
		return false;
	else
		return true;

/*	
	for(int i = 0; i < rows; i++){
			 
			if(agileMapUtils->SrcDist(i,0,l,b) < radius || agileMapUtils->SrcDist(i,cols,l,b) < radius) {
				return false;
			}
		}
	for(int j=0; j < cols; j++){
			if(agileMapUtils->SrcDist(0,j,l,b)<radius || agileMapUtils->SrcDist(rows,j,l,b)<radius) {
				return false;
			}
		}
	
	return true;
*/
 }

