


#ifndef _PLOT_CTS_2D_
#define _PLOT_CTS_2D_





void PlotCts2D(
	bool        shiftnorth_b,
	int         skysegmentation,
	const char* file,
	float       binsize,
	double      sourceradiusremove = 0,
	int         smoothing = 3,
	int         NCONREGSEARCH = 2,
	const char* outFile = "",
	bool        conregsearchType = 1,
	int         removeStep = 0,
	int         nstep = 1,
	int         bkgint = 20,
	const char* filesubtract = "",
	double      maxMapValue = -1);


void PlotCts2D_FINE(
	bool        shiftnorth_b,
	int         skysegmentation,
	const char* inputfile,
	float       binsize,
	double      sourceradiusremove = 0,
	int         smoothing = 3,
	int         NCONREGSEARCH = 2,
	const char* outFile = "",
	double      nearradious = 1,
	const char* expfile = "",
	double      minExp = 200,
	bool        conregsearchType = 1,
	int         removeStep = 0,
	int         nstep = 1,
	int         bkgint = 20,
	const char* filesubtract = "",
	double      maxMapValue = -1);




#endif

