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





ExpRatioEvaluator::ExpRatioEvaluator(	const char * _expPath,
					bool _isExpMapNormalized, 
					bool _createExpNormalizedMap,
					bool _createExpRatioMap, 
					double _minThreshold, 
					double _maxThreshold,
					double _squareSize
				    ) : agileMap(_expPath)
{	
	cout << "\n**ExpRatioEvaluator - started " << _isExpMapNormalized << " " << _createExpNormalizedMap << " " << _minThreshold << " " << _maxThreshold << " " << _squareSize << endl;
	
	isExpMapNormalized = _isExpMapNormalized;
	
	createExpNormalizedMap = _createExpNormalizedMap;

	createExpRatioMap = _createExpRatioMap;
	
	minThreshold = _minThreshold;
	
	maxThreshold = _maxThreshold;

	expPath=_expPath;
	
	cdelt2 = agileMap.GetYbin();	
 	squareSize = _squareSize/cdelt2;

	

	/*
		Reading from file .exp
	*/
		
	//image = MapConverter::mapPathToDoublePtr();
	//MapConverter::mapPathToDoublePtr(expPath);
	if(! convertFitsDataToMatrix() )
	{
		fprintf( stderr, "*ERROR: convertFitsDataToMatrix(): error reading fits file\n");
		exit (EXIT_FAILURE);
	}
	
	createAndWriteImages();
	
} 



ExpRatioEvaluator::ExpRatioEvaluator(
					AgileMap _agileMap, 
					bool _isExpMapNormalized, 
					bool _createExpNormalizedMap, 
					bool _createExpRatioMap, 
					double _minThreshold, 
					double _maxThreshold, 
					double _squareSize
					) : agileMap(_agileMap)
{

	cout << "\n**ExpRatioEvaluator - started " << _isExpMapNormalized << " " << _createExpNormalizedMap << " " << _minThreshold << " " << _maxThreshold << " " << _squareSize << endl;

	isExpMapNormalized = _isExpMapNormalized;
	
	createExpNormalizedMap = _createExpNormalizedMap;

	createExpRatioMap = _createExpRatioMap;
	
	minThreshold = _minThreshold;
	
	maxThreshold = _maxThreshold;

 



	/*
		Reading from AgileMap object 
	*/

 
	cdelt2 = agileMap.GetYbin();	
 	squareSize = _squareSize/cdelt2;
 	expPath = agileMap.GetFileName();
	rows = agileMap.Rows(); 
	cols = agileMap.Cols();



	image = new double*[rows];
 


	for (int i = 0; i < rows; ++i){
		image[i] = new double[cols];
 
	}
	
	double * data = agileMap.Buffer();

 
 
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols ; j++){
			image[cols-1-j][i] = data[i*cols+j];
		}
	}
	 

	createAndWriteImages();
	

	
}


 
void ExpRatioEvaluator::clearHeap(){

	cout << "*Clearing Heap.." << endl;
	// deleting image
	for (int i=0; i<cols; i++)
		delete [] image[i];
	delete [] image;


	// deleting normalized image (only if it has been created)
	if(! isExpMapNormalized){
		for (int i=0; i<cols; i++)
			delete [] normalizedImage[i];
		delete [] normalizedImage;

		for (int i=0; i<cols; i++)
    			delete [] normalizationFactorMatrix[i];
		delete [] normalizationFactorMatrix;
	}

	// deleting exp ratio image
	if(createExpRatioMap){
		for (int i=0; i<cols; i++)
			delete [] expRatioImage[i];
		delete [] expRatioImage;	
	}
	cout << "*Clearing Heap completed!" << endl;
	
}




void ExpRatioEvaluator::createAndWriteImages(){


	//if the exp map given in input is NOT already normalized 
	if(! isExpMapNormalized){
		normalizedImage = createNormalizedImage();

	}else{
		normalizedImage = image;
	}
		
	//expRatioNormalized image
	if(createExpNormalizedMap){
		writeMatrixDataInAgileMapFile("norm.exp", normalizedImage);
	}
	

	//expRatio image
	if(createExpRatioMap){
		expRatioImage = createExpRatioPixelMap();	
		writeMatrixDataInAgileMapFile("exp_norm.exp", expRatioImage);
	}
		

}

double ** ExpRatioEvaluator::getImage() {
	return image;
}

double ** ExpRatioEvaluator::getExpRatioMap(){
	return expRatioImage;
}

double ** ExpRatioEvaluator::getNormalizedMap(){
	return normalizedImage;
}

const char* ExpRatioEvaluator::getExpPath(){
	return expPath;
}
double ExpRatioEvaluator::getMinThreshold(){
	return minThreshold;
}
double ExpRatioEvaluator::getMaxThreshold(){
	return maxThreshold;
}
double ExpRatioEvaluator::getSquareSize(){
	cout << "ss: " << squareSize << endl;
	return squareSize;
}

bool ExpRatioEvaluator::convertFitsDataToMatrix()
{
	
	//CFITSIO
	fitsfile *fptr;    //FITS file pointer, defined in fitsio.h 
	int status = 0;   // CFITSIO status value MUST be initialized to zero! 
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
				printf("*ERROR: only 1D or 2D images are supported\n");
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

				// get memory for 1 row 
				pixels = (double *)malloc(naxes[0] * sizeof(double));

				if (pixels == NULL)
				{
					printf("*Memory allocation error - Press any key to exit\n");
					return false;
				}
				else
				{
					// loop over all the rows in the image, top to bottom 

					int col_index = 0;
					int row_index = 0;
					for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
					{
						if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status))  // read row of pixels 
							break;  // jump out of loop on error 

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
		cout << "*ERROR: Can't open fits file " << expPath << endl; 
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
	if(distSx < squareSize || distDx < squareSize || distUp < squareSize || distDown < squareSize)
		return false;
	else
		return true;
}




double ExpRatioEvaluator::computeExpRatioValues(double l, double b) 
{ 
	int x;
	int y;

	agileMap.GetRowCol(l,b,&x,&y);
	if(x < 0 || x > cols-1 || y < 0 || y > rows-1)
	{
		fprintf( stderr, "*ERROR: computeExpRatioValues(double l, double b):  Map .cts and Map .exp are not centered in the same galactic coordinate because input l and input b go out of the .exp map\n");
		output[0] =  -2;
		output[1] =  -2;
		output[2] =  -2;
		output[3] =  -2;
		
		return output[0];
	}
	return computeExpRatio(x, y);
}


double ExpRatioEvaluator::computeExpRatioValues(int i, int j, string type)
{	
	return computeExpRatio(i, j);
}


double ExpRatioEvaluator::computeExpRatio(int x, int y){
	
	// The output array  [ exp-ratio, nBad, nTot, greyLevelMean ]	
	double output[4];

 	int xmin, xmax, ymin,ymax;
	int npixel = 0;
	int nBad = 0;
	double totCount=0;
	double tmp = 0;
	double expRatio = 0;
	double greyLevelSum = 0;
	double mean = 0;
		
	xmin = x - squareSize;
	xmax = x + squareSize;
	ymin = y - squareSize;
	ymax = y + squareSize;


	if(isRectangleInside(x,y)) 
	{	
		for(int i = xmin; i <= xmax; i++) 
		{
			for(int j= ymin; j <= ymax; j++) 
			{	
				totCount+=1;
				//if(isExpMapNormalized) 
				tmp = (double)normalizedImage[i][j];
				//else
				//	tmp=(double)image[i][j];

				greyLevelSum+=tmp;
 				
				//cout << "minThreshold: " << minThreshold << "maxThreshold: " << maxThreshold << "tmp: " << tmp << endl;
				if(tmp < minThreshold || tmp >maxThreshold){
 					
					nBad+=1;
				}
					
			}
		}
		//cout << "nBad: " << nBad << "totCount: " << totCount << "expratio: " << (1-(nBad/totCount))*100 << endl;
 		//getchar();
		output[0] = (1-(nBad/totCount))*100;
		output[1] = nBad;
		output[2] = totCount;
		output[3] = greyLevelSum/totCount;
		

		return output[0];

	}else 
	{
		//fprintf( stderr, "Rectangle is not completely inside!\n");
		output[0] =  -1;
		output[1] =  -1;
		output[2] =  -1;
		output[3] =  -1;

		return output[0];

	}
}



double ** ExpRatioEvaluator::createExpRatioPixelMap(){


		cout << "*Creating exp ratio image.."<<endl;
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
 
				expRatioMap[i][j] = computeExpRatioValues(i,j,"PIXEL");
								
			}
			
		}
		
		cout << "*Creating exp ratio image completed!"<<endl;
		return expRatioMap;
}




double ExpRatioEvaluator::Alikesinaa(double input){
	
		if (input == 0)
			 
			return(1.0);
				
		else
			return(sin(input)/input);
}

	
	

// OPTIMIZED CODE!
double ** ExpRatioEvaluator::computeSpatialNormalizationFactorMatrix(){
			
			//cout << "*Computing normalization factors matrix.." << endl;
			/*
					D D D X X X
					D D D X X X
					D D D X X X
					V V V O O O
					V V V O O O
					V V V O O O
			*/		
			

			// Initializes a distance matrix.
			int dMrows = rows/2;
			int dMcols =  cols/2;
			double center_l = agileMap.GetMapCenterL();
			double center_b = agileMap.GetMapCenterB();
			double fctr4Normalization = 0.0003046174197867085688996857673060958405 * cdelt2 * cdelt2;
			
			 

			double ** normalizationFactorMatrix = new double*[rows];
			for(int i = 0; i < rows; ++i) {
				normalizationFactorMatrix[i] = new double[cols];
			}
	
		/*		
			Computes distance from double ** image using SrcDist( i, j, center_l, center_b ) to store the distances from every pixel of
			the first quadrant of the matrix (i,j) to the center (center_l, center_b).

			Initialize and compute the normalizationFactorMatrix, iterating ONLY ON THE FIRST QUADRANT OF THE MATRIX, WRITING 4 PIXEL (the symmetric 				ones) AT TIME, with this formula:

			if the distance is not 0 
			0.0003046174197867085688996857673060958405 * xbin * ybin * Alikesinaa(0.0174532925199432954743716805978692718782 * distanceMatrix[i][j]);
			else
			0.0003046174197867085688996857673060958405 * xbin * ybin

		*/	
		
				
		for(int i=0; i<dMrows; ++i) {
			for(int j=0; j<dMcols; ++j) {
				
				double distance =agileMap.SrcDist(i,j,center_l,center_b);

 
 

				if(distance!=0){
				
					normalizationFactorMatrix[i][j] =  fctr4Normalization * Alikesinaa(0.0174532925199432954743716805978692718782 * distance);	 
					
					normalizationFactorMatrix[i][cols -1 - j] =fctr4Normalization * Alikesinaa(0.0174532925199432954743716805978692718782 * distance);
					
					normalizationFactorMatrix[rows -1 - i][j] = fctr4Normalization * Alikesinaa(0.0174532925199432954743716805978692718782 * distance);
					
					normalizationFactorMatrix[rows -1 - i][cols -1 - j] = fctr4Normalization * Alikesinaa(0.0174532925199432954743716805978692718782 * distance);
					
				}
				
				else {
					normalizationFactorMatrix[i][j] = fctr4Normalization;
					normalizationFactorMatrix[i][cols - j] = fctr4Normalization;
					normalizationFactorMatrix[rows - i][j] = fctr4Normalization;
					normalizationFactorMatrix[rows - i][cols - j] = fctr4Normalization;
				}
			}
		}

		//cout << "*Computing normalization factors matrix completed!" << endl;
		return normalizationFactorMatrix;
				
}



double ** ExpRatioEvaluator::createNormalizedImage(){

	//cout << "*Normalizing image... " << endl;
	
// Computes time normalization factor
	double timeFactor = agileMap.GetTstop()  -  agileMap.GetTstart();
 
	// Computes spatial normalization factor matrix
	double ** normalizationFactorMatrix = computeSpatialNormalizationFactorMatrix();

	
	// Computes normalizedImage!! 
	double ** normalizedImage = new double*[rows];
	for(int i = 0; i < rows; ++i) {
		normalizedImage[i] = new double[cols];
	}
			
	for(int i = 0 ; i < rows; ++i) {
		for(int j = 0 ; j < cols; ++ j) {
			normalizedImage[i][j] = image[i][j]/ (  timeFactor * normalizationFactorMatrix[i][j] );

		}
	}
			
	//cout << "*Normalization completed!" << endl;
	return normalizedImage;

}
 
 
int ExpRatioEvaluator::writeMatrixDataInAgileMapFile(const char * appendToFilename, double ** matrixData){  
	
	cout << "*Writing data " << appendToFilename <<" in fits file.."<<endl;

	/// Computes new filename
	const char * newFileNameC;

	string cleanedFilename = computeNewFileName(appendToFilename); // TO FIX 
 	//cout << "cleanedFilename: " << cleanedFilename << endl;


	newFileNameC = cleanedFilename.c_str();
	//cout << "new file name c_str= " << newFileNameC << endl;

        remove(newFileNameC);

 	FitsFile f;

	if (!f.Create(newFileNameC)) {
		cerr << "*ERROR " << f.Status() << " creating " << newFileNameC << endl;
		return f.Status();
	}
	else
		cout << "*Created "<<newFileNameC << endl;	

	
	int bitpix = DOUBLE_IMG;
	int naxis = 2;
	long naxes[2] = { rows, cols };		 
	long fpixel[2] = { 1, 1 };

	fits_create_img(f, bitpix, naxis, naxes, f);
	
	double * rowData =  new double[rows*cols];
	
	int status;
	
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols ; j++){
			rowData[i*cols+j] = matrixData[rows-i-1][j];
		}
	}	
	 	
	
	fits_write_pix(f, TDOUBLE, fpixel, rows*cols, rowData, f);

	f.WriteKey("CTYPE1", "GLON-ARC");
	f.WriteKey("CTYPE2", "GLAT-ARC");
	f.WriteKey("CRPIX1", agileMap.GetX0());
	f.WriteKey("CRVAL1", agileMap.GetMapCenterL());
	f.WriteKey("CDELT1", agileMap.GetXbin());
	f.WriteKey("CRPIX2", agileMap.GetY0());
	f.WriteKey("CRVAL2", agileMap.GetMapCenterB());
	f.WriteKey("CDELT2", agileMap.GetYbin());
	f.WriteKey("LONPOLE", agileMap.GetLonpole());
	f.WriteKey("MINENG", agileMap.GetEmin());
	f.WriteKey("MAXENG", agileMap.GetEmax());
	f.WriteKey("INDEX", agileMap.GetMapIndex());
	f.WriteKey("SC-Z-LII", agileMap.GetLpoint());
	f.WriteKey("SC-Z-BII", agileMap.GetBpoint());
	//f.WriteKey("SC-LONPL", m_gp);

	const char * m_dateObs = agileMap.GetStartDate();
	const char * m_dateEnd = agileMap.GetEndDate();

	if (m_dateObs[0])
		f.WriteKey("DATE-OBS", m_dateObs);
	if (m_dateEnd[0])
		f.WriteKey("DATE-END", m_dateEnd);

	f.WriteKey("TSTART", agileMap.GetTstart());
	f.WriteKey("TSTOP", agileMap.GetTstop());
	f.WriteKey("FOVMIN", agileMap.GetFovMin());
	f.WriteKey("FOV", agileMap.GetFovMax());
	f.WriteKey("ALBEDO", agileMap.GetAlbedo());
	f.WriteKey("PHASECOD", agileMap.GetPhaseCode());

	double m_step = agileMap.GetStep();
	if (m_step)
		f.WriteKey("STEP", m_step);

	const char * m_skyL = agileMap.GetSkyL();
	const char * m_skyH = agileMap.GetSkyH();

	if (m_skyL[0])
		f.WriteKey("SKYL", m_skyL);
	if (m_skyH[0])
		f.WriteKey("SKYH", m_skyH);

	if (f.Status())
		cerr << "*ERROR " << f.Status() << " writing to " << newFileNameC << endl;
	else
		cout << "*Writing data to" << newFileNameC <<" completed!"<<endl;
 

	return f.Status();


	

}

string ExpRatioEvaluator::computeNewFileName(const char * appendToFilename){
	
	cout << "*Computing fits filename.."<<endl;

	//cout << "append: " << appendToFilename << endl;
	
  	const char * imageName = agileMap.GetFileName(); // e.g.    MAP1000s.exp
	
	//cout << "imageName_cstr: " << imageName << endl;

	string imageName_string(imageName);
	

	//cout << "imageName_str: " << imageName << endl;

	string newFileName = "";

	// tolgo tutto quello che c'è prima di tutti gli slash	
	size_t foundPatternSlash = imageName_string.find("/");
	while(foundPatternSlash != string::npos)
	{ 
		imageName_string = imageName_string.substr(foundPatternSlash+1);
		//cout << "imageName senza /: " << imageName_string << endl;
		foundPatternSlash = imageName_string.find("/");
	}



	
	
	size_t foundPatternExp = imageName_string.find(".exp");
	size_t foundPatternCts = imageName_string.find(".cts");

	if(foundPatternExp != string::npos)
		newFileName = imageName_string.substr(0,foundPatternExp);
	else if(foundPatternCts != string::npos)
		newFileName = imageName_string.substr(0,foundPatternCts);
	else
		newFileName = imageName_string;

	//cout << "new file name = " << newFileName << endl;

	string appendToFilenameString(appendToFilename);

	// convert double to string
	ostringstream minThresholdStringStream;
	minThresholdStringStream << minThreshold;
	string min_str = minThresholdStringStream.str();
	size_t foundPatternDotMin = min_str.find(".");
	min_str = min_str.substr(0,foundPatternDotMin);		
	

	// convert double to string
	ostringstream maxThresholdStringStream;
	maxThresholdStringStream << maxThreshold;
	string max_str = maxThresholdStringStream.str();
	size_t foundPatternDotMax = max_str.find(".");
	max_str = max_str.substr(0,foundPatternDotMax);		

	//string max_str = to_string(maxThreshold);
	//size_t foundPatternDotMax = max_str.find(".");
	//max_str = max_str.substr(0,foundPatternDotMax);		
	
	
	ostringstream sqrSizeStringStream;
	sqrSizeStringStream<<squareSize*cdelt2;
	string sqrSize_str = sqrSizeStringStream.str();
    	
    
	   

   	newFileName +="_"+min_str+"_"+max_str+"_"+sqrSize_str+"_"+appendToFilenameString+".gz";		
   
	//cout << "newFileName: " << newFileName << endl;  

	//replace( newFileName.begin(), newFileName.end()-7, '.', '_');
	// cout << "NewFileName: " << newFileName << endl;


	cout << "*Computing fits filename completed!"<<endl;

	return newFileName;
}


