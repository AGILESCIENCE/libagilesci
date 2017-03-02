/*
 * Copyright (c) 2005-2016
 *     Andrew Chen, Alberto Pellizzoni, Alessio Trois (IASF-Milano),
 *     Andrea Bulgarelli, Andrea Zoli (IASF-Bologna),
 *     Tomaso Contessi (Nuove Idee sas)
 *
 * Any information contained in this software
 * is property of the AGILE TEAM and is strictly
 * private and confidential.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <prj.h>
#include <wcstrig.h>
#include <sph.h>
#include <CalibUtils.h>
#include <FitsUtils.h>
#include <MathUtils.h>

#include "Selection.h"
#include "Eval.h"

//#define DEBUG 1

using std::cout;
using std::cerr;
using std::endl;

namespace eval
{

const double obtlimit = 104407200.0;
enum ProjectionType {ARC, AIT};

struct MapspecEntry
{
    double fovradmin;
    double fovradmax;
    double emin;
    double emax;
    double index;
};

class Mapspec : public std::vector<MapspecEntry>
{
public:
    Mapspec() {}
    bool Store(const char* maplist);
    void Print();
};

bool Mapspec::Store(const char* maplist)
{
    bool stillgoing = true;
    bool succeeded = false;
    MapspecEntry mapspec;
    std::fstream mapstream(maplist);
    while (stillgoing) {
        mapstream >> mapspec.fovradmin >> mapspec.fovradmax >> mapspec.emin >> mapspec.emax >> mapspec.index;
        if (mapstream.good()) {
            push_back(mapspec);
            succeeded = true;
        } else
            stillgoing = false;
    }
    return succeeded;
}

void Mapspec::Print()
{
    cout << "Energy ranges (FOVmin [deg], FOVmax [deg], Emin [MeV], Emax [MeV], spectral index) :" << endl;
    for (iterator i = begin(); i != end(); i++)
        cout << std::setprecision(2) << "     "
             << i->fovradmin << "   "
             << i->fovradmax << "     "
             << i->emin << "   "
             << i->emax << "   "
             << i->index << endl;
}


double Alikesinaa(double input)
{
    if (input == 0) return(1.0); else return(sin(input)/input);
}

double AG_expmapgen_area(double xbin, double ybin, double theta, ProjectionType proj)
{
    if (proj==ARC)
        return 0.0003046174197867085688996857673060958405 * xbin * ybin * Alikesinaa(0.0174532925199432954743716805978692718782 * theta);
    else
        return 0.0003046174197867085688996857673060958405 * xbin * ybin;
}

void ReadFitsCol(fitsfile* file, double* array, const char* colName, long rowOffs, long rowCount, int* status)
{
    int colNum;
    fits_get_colnum(file, 1, const_cast<char*>(colName), &colNum, status);
    fits_read_col(file, TDOUBLE, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, status);
}

void ReadFitsCol(fitsfile* file, short* array, const char* colName, long rowOffs, long rowCount, int* status)
{
    int colNum;
    fits_get_colnum(file, 1, const_cast<char*>(colName), &colNum, status);
    fits_read_col(file, TSHORT, colNum, rowOffs+1, 1, rowCount, NULL, array, NULL, status);
}

bool inmap(int i, int ii, int mxdim)
{
    return (i < mxdim) && (i >= 0) &&
           (ii < mxdim) && (ii >= 0);
}

bool LoadTimeList(const char* timelist, Intervals& intervals, double& tmin, double& tmax)
{
    if (strcmp(timelist, "None")) {
        intervals = ReadIntervals(timelist);
        int count = intervals.Count();
        if (!count)
            return false;
        tmin = intervals.Min();
        tmax = intervals.Max();
        cout << "Time intervals:" << endl;
        for (int i=0; i<count; i++)
            cout << intervals[i].Start() << " " << intervals[i].Stop() << endl;
    }
    else {
        Interval intv(tmin, tmax);
        intervals.Add(intv);
    }
    if (tmin <= obtlimit) {
        cerr << "ERROR: initial time must be greater than " << obtlimit << endl;
        return false;
    }

    return true;
}

bool AlbTest(double ra, double dec, double earth_ra, double earth_dec, double albrad)
{
    return SphDistDeg(ra, dec, earth_ra, earth_dec) > albrad;
}

bool FovTest(Mapspec &maps, long i, double theta)
{
    return theta < maps[i].fovradmax && theta >= maps[i].fovradmin;
}

void addexpval(long i, long ii, long mxdim, const vector<int> &aitstatus,
        const vector<double>& lng, const vector<double>& lat, double lp, double bp,
        Mapspec &maps, double learth, double bearth, double albrad, double time, vector<double> &A,
        AeffGridAverage *raeffArr, const vector<double>& area)
{
    long element = i * mxdim + ii;
    if (aitstatus[element] == 0 && AlbTest(lng[element], lat[element], learth, bearth, albrad)) {
        double theta = SphDistDeg(lng[element], lat[element], lp, bp);
        double phi = 0.0;
        long nmaps = maps.size();
        for (long m = 0; m < nmaps; m++) {
            if (FovTest(maps, m, theta)) {
                A[m*mxdim*mxdim+element] += 1e-3*time*(raeffArr[m].AvgVal(theta, phi))*area[element];
            }
        }
    }
}

bool LoadProjection(const char *projstr, ProjectionType &proj)
{
    if ((strcmp(projstr, "ARC") == 0) || (strcmp(projstr, "arc") == 0) ||
        (strcmp(projstr,"GLON-ARC") == 0) || (strcmp(projstr,"GLAT-ARC") == 0))
        proj = ARC;
    else if ((strcmp(projstr, "AIT") == 0) || (strcmp(projstr, "ait") == 0) ||
             (strcmp(projstr, "AITOFF") == 0) || (strcmp(projstr,"aitoff") == 0) ||
             (strcmp(projstr, "GLON-AIT") == 0) || (strcmp(projstr, "GLAT-AIT") == 0))
        proj = AIT;
    else {
        cerr << endl << "ERROR: Invalid Projection" << endl;
        return false;
    }

    return true;
}

bool YTolTest(double ra_y, double dec_y, double ra_y0, double dec_y0, double y_tol)
{
    return SphDistDeg(ra_y, dec_y, ra_y0, dec_y0) > y_tol;
}

bool EarthTolTest(double earth_ra, double earth_dec, double earth_ra0, double earth_dec0, double earth_tol)
{
    return SphDistDeg(earth_ra, earth_dec, earth_ra0, earth_dec0) > earth_tol;
}

const char *AfterSlashPtr(const char *str)
{
    const char *slash = strrchr(str, '/');
    return slash ? slash+1 : str;
}

int WriteTime(fitsfile* input, double t1, double t2)
{
    const double mjdrefi = 53005.;
    const double mjdreff_OLD = 0.00075444444444444444;
    const double mjdreff_2009 = 0.00076601852;
    const double mjdreff_2012 = 0.000777593;

    int status = 0;
    int hdnum = 0;
    fits_get_num_hdus(input, &hdnum, &status);

    int i;
    fits_get_hdu_num(input, &i);

    char str[100] = "TT";
    float timezero = 0.0;
    char secstr[] = "s";

    fits_update_key(input, TDOUBLE, "TSTART", &t1,  "[OBT]first event time", &status);
    fits_update_key(input, TDOUBLE, "TSTOP", &t2,  "[OBT]last event time", &status);
    fits_update_key(input, TSTRING, "TIMESYS",  str, NULL, &status);
    fits_update_key(input, TSTRING, "TIMEUNIT",  secstr, NULL, &status);

    double mjdreff = mjdreff_OLD;
    if (t1 >= 268185600.0)
        mjdreff = mjdreff_2012;
    else if (t1 >= 157852800.0)
        mjdreff = mjdreff_2009;

    int oldstatus = status;
    time_t tt;
    time(&tt);
    tm* tstr;
    tstr = gmtime(&tt);

    tstr->tm_year = 107;
    tstr->tm_mon = 0;
    tstr->tm_mday = 1;
    tstr->tm_hour = 0;
    tstr->tm_min = 0;
    tstr->tm_sec = 0;
    tstr->tm_isdst = 0;
    tt = timegm(tstr);
    double tinc = difftime(tt+1,tt);

    tt += time_t(((mjdrefi - 54101 + mjdreff) * 86400 + timezero + t1) / tinc);
    tstr = gmtime(&tt);
    double inttstart;
    double sec = tstr->tm_sec + modf(t1, &inttstart);
    char datestr[FLEN_KEYWORD];
    fits_time2str(tstr->tm_year + 1900, tstr->tm_mon + 1, tstr->tm_mday, tstr->tm_hour, tstr->tm_min, sec, 0, datestr, &status);
    fits_update_key(input, TSTRING, "DATE-OBS",  datestr, "start date and time of the observation(TT)", &status);

    tt += (time_t)((t2 - inttstart) / tinc);
    tstr = gmtime(&tt);
    sec = tstr->tm_sec + modf(t2, &inttstart);
    char datestr1[FLEN_KEYWORD];
    fits_time2str(tstr->tm_year + 1900, tstr->tm_mon + 1, tstr->tm_mday, tstr->tm_hour, tstr->tm_min, sec, 0, datestr1, &status);
    fits_update_key(input, TSTRING, "DATE-END",  datestr1, "end date and time of the observation(TT)", &status);

    if (i == hdnum) {
        fits_movabs_hdu(input, 1, NULL, &status);
        fits_update_key(input, TDOUBLE, "TSTART", &t1,  "[OBT]first event time", &status);
        fits_update_key(input, TDOUBLE, "TSTOP", &t2,  "[OBT]last event time", &status);
        fits_update_key(input, TSTRING, "TIMESYS",  str, NULL, &status);
        fits_update_key(input, TSTRING, "TIMEUNIT",  secstr, NULL, &status);
        fits_update_key(input, TFLOAT, "TIMEZERO",  &timezero, str, &status);
        fits_update_key(input, TSTRING, "DATE-OBS",  datestr, "start date and time of the observation(TT)", &status);
        fits_update_key(input, TSTRING, "DATE-END",  datestr1, "end date and time of the observation(TT)", &status);
        fits_movabs_hdu(input, i, NULL, &status);
        }

    for (int i=1; i <= hdnum; i++) {
        fits_movabs_hdu(input, i, NULL, &status);
        oldstatus = status;
        fits_delete_key(input, "MJDREF", &status);
        if (status == 202)
            status = oldstatus;
        double mjdrefiVal = mjdrefi;
        fits_update_key(input, TDOUBLE, "MJDREFI",  &mjdrefiVal, NULL, &status);
        fits_update_key(input, TDOUBLE, "MJDREFF",  &mjdreff, NULL, &status);
    }

    return status;
}

int EvalExposure(const char *outfile, const char *sarFileName,
                 const char *edpFileName, const char *maplist,
                 const char *projection, double mdim, double mres,
                 double la, double ba, double lonpole, double albrad,
                 double y_tol, double roll_tol, double earth_tol,
                 int phasecode, double binstep, int timestep,
                 double index, double tmin, double tmax, double emin,
                 double emax, double fovradmin, double fovradmax,
                 const char *selectionFilename, const char *templateFilename,
                 Intervals &intervals, vector< vector<double> > &exposures, bool saveMaps)
{
    int status = 0;

    Mapspec maps;
    if (strcmp(maplist, "None")) {
        maps.Store(maplist);
        maps.Print();
    }
    else {
        MapspecEntry mapspec;
        mapspec.fovradmin = fovradmin;
        mapspec.fovradmax = fovradmax;
        mapspec.emin = emin;
        mapspec.emax = emax;
        mapspec.index = index;
        maps.push_back(mapspec);
    }
    long nmaps = maps.size();

    long mxdim = long(mdim / mres + 0.1); // dimension (in pixels) of the map
#ifdef DEBUG
    cout << "mdim: " << mdim << " mres: " << mres << " mxdim: " << mxdim << " nmaps: " << nmaps << endl;
#endif
    long npixels = mxdim * mxdim;
    exposures.resize(intervals.Count());
    for (int i=0; i < intervals.Count(); i++) {
        exposures[i].resize(nmaps * npixels);
        for (int j=0; j < nmaps * npixels; j++)
            exposures[i][j] = 0;
    }

    struct prjprm *prj = new(prjprm);
    prjini(prj);
    ProjectionType proj;
    if (!LoadProjection(projection, proj)) {
        cerr << "Error loading projection '" << projection << "'" << endl;
        return false;
    }

    double lng0, lat0, x, y;
    long i, ii;
    vector<double> lng(npixels);
    vector<double> lat(npixels);
    vector<double> area(npixels);
    vector<int> aitstatus(npixels);
    switch (proj) {
    case ARC:
        double theta, phi;
        double eul[5];
        for (i = 0; i <= mxdim-1; i++) {
            x = -(mdim/2)+(mres*(i+0.5));
            for (ii = 0; ii <= mxdim-1; ii++) {
                aitstatus[i*mxdim+ii] = 0;
                y = -(mdim/2)+(mres*(ii+0.5));
                theta = 90-sqrt(x*x+y*y);
                phi = atan2d(-y, -x);
                eul[0] = la;
                eul[1] = 90-ba;
                eul[2] = lonpole;
                eul[3] = cosd(eul[1]);
                eul[4] = sind(eul[1]);
                sphx2s(eul, 1, 1, 0, 0, &phi, &theta, &lng0, &lat0);
                lng[i*mxdim+ii] = lng0;
                lat[i*mxdim+ii] = lat0;
                area[i*mxdim+ii] = AG_expmapgen_area(mres, mres, 90-theta, proj);
#ifdef DEBUG
                cout << "lng[" << i << "][" << ii << "]=" << lng[i*mxdim+ii] << endl;
                cout << "lat[" << i << "][" << ii << "]=" << lat[i*mxdim+ii] << endl;
                cout << "area[" << i << "][" << ii << "]=" << area[i*mxdim+ii] << endl;
#endif
            }
        }
        break;
    case AIT:
        int ait0;
        prj->r0 = 180/M_PI;
        prj->phi0 = la;
        prj->theta0 = ba;
        aitset(prj);
        for (i = 0; i <= mxdim-1; i++) {
            y = -(mdim/2)+(mres*(i+0.5));
            for (ii = 0; ii <= mxdim-1; ii++) {
                x = (mdim/2)-(mres*(ii+0.5));
                aitx2s(prj, 1, 1, 0, 0, &x, &y, &lng0, &lat0, &ait0);
                lng[i*mxdim+ii] = lng0;
                lat[i*mxdim+ii] = lat0;
                aitstatus[i*mxdim+ii] = ait0;
                area[i*mxdim+ii] = AG_expmapgen_area(mres, mres, 0, proj);
            }
        }
        break;
    }

    long n = 0;
    double time = 0;

    int hdutype = 0;
    int  find = 0;

    const int size = 65536;
    double ra_y[size];
    double dec_y[size];
    double livetime[size];
    double earth_ra[size];
    double earth_dec[size];
    short  phase[size];

    AeffGridAverage *raeffArr = new AeffGridAverage[nmaps];

    bool hasEdp = strcmp(edpFileName, "None");
    for (long m = 0; m < nmaps; m++) {
        int status = raeffArr[m].Read(sarFileName, maps[m].emin, maps[m].emax, maps[m].index);
        if (status) {
            cerr << "Error reading " << sarFileName << endl;
            delete []raeffArr;
            return status;
        }
        if (hasEdp)
            status = raeffArr[m].LoadEdp(edpFileName);
        if (status) {
            cerr << "Error reading " << edpFileName << endl;
            delete []raeffArr;
            return status;
        }
    }

    fitsfile* selectionFits;
    if (fits_open_file(&selectionFits, selectionFilename, READONLY, &status)) {
        cerr << "ERROR opening selection file " << selectionFilename << endl;
        return -1;
    }
    fits_movabs_hdu(selectionFits, 2, &hdutype, &status);

    fitsfile* templateFits;
    if (fits_open_file(&templateFits, templateFilename, READWRITE, &status)) {
        cerr << "ERROR opening template file " << templateFilename << endl;
        return -1;
    }
    fits_movabs_hdu(templateFits, 2, &hdutype, &status);
    long oldnrows;
    fits_get_num_rows(templateFits, &oldnrows, &status);
    fits_delete_rows(templateFits, 1, oldnrows, &status);

#ifdef DEBUG
    cout << "Evaluating exposure.." << endl;
#endif

    char tempFilename[FLEN_FILENAME];
    tmpnam(tempFilename);
    double lp = 0., bp = 0.;
    double learth, bearth;
    double lp0 = 0., bp0 = 0., gp0 = 0.;
    for (int intvIndex = 0; intvIndex < intervals.Count(); intvIndex++) {
#ifdef DEBUG
        cout << "Interval #" << intvIndex << endl;
#endif
        vector<double> &A = exposures[intvIndex];
        Intervals sIntv;
        sIntv.Add(intervals[intvIndex]);
        string selExpr = selection::TimesExprString(sIntv);
#ifdef DEBUG
        cout << "selExpr: " << selExpr << endl;
#endif
        fits_select_rows(selectionFits, templateFits, (char*)selExpr.c_str(), &status);
#ifdef DEBUG
        cout << "Rows from " << tempFilename << " selected" << endl;
#endif
        find++;

        long allnrows;
        fits_get_num_rows(templateFits, &allnrows, &status);
#ifdef DEBUG
        cout << "all nrows: " << allnrows << endl;
#endif

        long rowblockincrement = size;
        fits_get_rowsize(templateFits, &rowblockincrement, &status);
#ifdef DEBUG
        cout << "Row block size = " << rowblockincrement << endl;
#endif
        long* change = new long[rowblockincrement];

        for (long rowblockzero=0; rowblockzero < allnrows; rowblockzero += rowblockincrement) {
            long nrows = (rowblockzero + rowblockincrement < allnrows) ? rowblockincrement : (allnrows - rowblockzero);

#ifdef DEBUG
            cout << "Reading " << nrows << " rows" << endl;
#endif
            ReadFitsCol(templateFits, ra_y, "ATTITUDE_RA_Y", rowblockzero, nrows, &status);
            ReadFitsCol(templateFits, dec_y, "ATTITUDE_DEC_Y", rowblockzero, nrows, &status);
            ReadFitsCol(templateFits, livetime, "LIVETIME", rowblockzero, nrows, &status);
            ReadFitsCol(templateFits, phase, "PHASE", rowblockzero, nrows, &status);
            ReadFitsCol(templateFits, earth_ra, "EARTH_RA", rowblockzero, nrows, &status);
            ReadFitsCol(templateFits, earth_dec, "EARTH_DEC", rowblockzero, nrows, &status);

            long count = 0;
            double earth_ra0 = earth_ra[0], earth_dec0 = earth_dec[0];
            double ra_y0 = ra_y[0], dec_y0 = dec_y[0];

            for (long k = 1; k<nrows; k++) {
                if ((phase[k-1] != phase[k]) || YTolTest(ra_y[k-1], dec_y[k-1], ra_y0, dec_y0, y_tol) || EarthTolTest(earth_ra[k-1], earth_dec[k-1], earth_ra0, earth_dec0, earth_tol)) {
                    change[count++] = k;
                    earth_ra0 = earth_ra[k-1];
                    earth_dec0 = earth_dec[k-1];
                    ra_y0 = ra_y[k-1];
                    dec_y0 = dec_y[k-1];
                }
            }
            if (count == 0)
                change[0] = nrows;
#ifdef DEBUG
            cout << endl << rowblockzero << " Count = " << count << endl;
#endif

            long lowrow = 0, highrow = 0;
            for (long k = 0; k<=count; k++) {
                time = 0;
                if (k == 0) {
                    lowrow = 0;
                    highrow = change[0];
                }
                else if (k == count) {
                    lowrow = change[k-1];
                    highrow = nrows;
                }
                else {
                    lowrow = change[k-1];
                    highrow = change[k];
                }

                for (n = lowrow; n<highrow; n++)
                    time += (livetime[n] * timestep);

                //euler(ra_y[lowrow], dec_y[lowrow], psi[lowrow], &lp, &bp, gp+lowrow);

                /// eulerold(ra_y[lowrow], dec_y[lowrow], &lp, &bp, 1);
                Euler(ra_y[lowrow], dec_y[lowrow], &lp, &bp, 1);

                /// eulerold(earth_ra[lowrow], earth_dec[lowrow], &learth, &bearth, 1);
                Euler(earth_ra[lowrow], earth_dec[lowrow], &learth, &bearth, 1);

                if (k == count && (rowblockzero + nrows) >= allnrows) {
                    lp0 = lp;
                    bp0 = bp;
                    if (isnan(lp0)) lp0 = 0.;
                    if (isnan(bp0)) bp0 = 0.;
                    //gp0=gp[lowrow]*R2D;
                }

                // reportfile << "LP = " << lp << " BP = " << bp << " LonGPole= " << R2D*gp[lowrow]<< endl;
                // reportfile << "Earth (theta, phi) = " << learth << ", " << bearth << endl;
                // cout << binstep << " " << mxdim << endl;
                /// Rotation satRot(q4[lowrow], q1[lowrow], q2[lowrow], q3[lowrow]);

                if (mxdim == 1) {
                    addexpval(0, 0, 1, aitstatus, lng, lat, lp, bp, maps, learth, bearth, albrad, time, A, raeffArr, area);
                }
                else {
                    for (i = 0; i <= mxdim-2; i += binstep) {
                        for (ii = 0; ii <= mxdim-2; ii += binstep)
                            addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, maps, learth, bearth, albrad, time, A, raeffArr, area);
                        ii = mxdim - 1;
                        addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, maps, learth, bearth, albrad, time, A, raeffArr, area);
                    }
                    i = mxdim - 1;
                    for (ii = 0; ii <= mxdim-2; ii += binstep)
                        addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, maps, learth, bearth, albrad, time, A, raeffArr, area);
                    ii = mxdim - 1;
                    addexpval(i, ii, mxdim, aitstatus, lng, lat, lp, bp, maps, learth, bearth, albrad, time, A, raeffArr, area);
                }
            }
        }

#ifdef DEBUG
        for(unsigned int m=0; m<maps.size(); m++) {
            for (i = 0; i <= mxdim-2; i+= binstep) {
                for (ii = 0; ii <= mxdim-2; ii+= binstep)
                    cout << A[m*mxdim*mxdim + i*mxdim+ii] << " ";
                cout << endl;
            }
            cout << endl;
        }
#endif
        delete []change;
        if (allnrows > 0)
            fits_delete_rows(templateFits, 1, allnrows, &status);
        if (status) {
            delete []raeffArr;
            return status;
        }
    }
    delete []raeffArr;
#ifdef DEBUG
    cout << "Ending evaluation" << endl;
#endif

    if (find == 0)
        return 1005;

    for (int intvIndex=0; intvIndex<intervals.Count(); intvIndex++) {
        Interval intv = intervals[intvIndex];
        vector<double> &A = exposures[intvIndex];

        long nrows = mxdim;
        long ncols = mxdim;
        long outpixel[2];
        double pixel = 0.0;
        // AB
        // 1 3
        // 4 2
        double pixel2 = 0.0;
        double pixel3 = 0.0;
        double pixel4 = 0.0;

        if(binstep > 1) {
            cout << "Beginning interpolation.." << endl;
            for (long long m = 0; m < nmaps; m++) {
                int step0 = binstep;
                bool end0 = true;
                for (outpixel[0]=1; end0; outpixel[0]+=binstep) {
                    bool end1 = true;
                    int step1 = binstep;
                    if (outpixel[0] + binstep > nrows) {
                        step0 = nrows - outpixel[0];
                        end0 = false;
                    }
                    for(outpixel[1]=1; end1; outpixel[1]+=binstep) {
                        if (outpixel[1] + binstep > ncols) {
                            step1 = ncols - outpixel[1];
                            end1 = false;
                        }
                        //cout << outpixel[1]-1 << "," << outpixel[0]-1 << ":" << end1 << "," << end0 << ":" << step1 << "," << step0 << ";";
                        pixel = A[m*npixels+(outpixel[1]-1)*mxdim+(outpixel[0]-1)];
                        if(pixel != 0.0 && step0 > 0 && step1 > 0) {
                            pixel2 = A[m*npixels+(outpixel[1]-1 + step1)*mxdim+(outpixel[0]-1 + step0)];
                            pixel3 = A[m*npixels+(outpixel[1]-1 + step1)*mxdim+(outpixel[0]-1)];
                            //si lavora sul primo triangolo
                            long y1 = outpixel[0];
                            long x1 = outpixel[1];
                            double z1 = pixel;
                            long y2 = outpixel[0] + step0;
                            long x2 = outpixel[1] + step1;
                            double z2 = pixel2;
                            long y3 = outpixel[0];
                            long x3 = outpixel[1] + step1;
                            double z3 = pixel3;
                            double a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
                            double b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
                            long c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
                            //cout << "1 " << x1 << " " << y1 << " " << z1 << endl;
                            //cout << "2 " << x2 << " " << y2 << " " << z2 << endl;
                            //cout << "3 " << x3 << " " << y3 << " " << z3 << endl;
                            for(int i=y1; i<=y2; i++) {
                                for(int j=x1; j<=x2; j++) {
                                    if(i-y1<=j-x1) {
                                        double p = (c * z1 - a * (j - x1) - b * (i - y1)) / c;
                                        long outpixel2[2];
                                        long x, y;
                                        outpixel2[0] = y = i;
                                        outpixel2[1] = x = j;
                                        if (x < 1 || x > mxdim || y < 1 || y > mxdim )
                                            cout << "x " << x << " " << y << " " << p << endl;
                                        A[m*npixels+(outpixel2[1]-1)*mxdim+(outpixel2[0]-1)] = p;
                                    }
                                }
                            }

                            pixel4 = A[m*npixels+(outpixel[1] - 1)*mxdim+(outpixel[0]-1 + step0)];
                            //cout << pixel << "," <<  pixel2 << "," <<  pixel3 << "," <<  pixel4 << "; ";
                            //si lavora sul secondo triangolo
                            y3 = outpixel[0] + step0;
                            x3 = outpixel[1];
                            z3 = pixel4;
                            a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
                            b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1);
                            c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
                            //cout << "4 " << x3 << " " << y3 << " " << z3 << endl;
                            for(int i=y1; i<=y2; i++) {
                                for(int j=x1; j<=x2; j++) {
                                    if(i-y1>j-x1) {
                                        double p = (c * z1 - a * (j - x1) - b * (i - y1)) / c;
                                        long outpixel2[2];
                                        long x, y;
                                        outpixel2[0] = y = i;
                                        outpixel2[1] = x = j;
                                        if (x < 1 || x > mxdim || y < 1 || y > mxdim )
                                            cout << "x " << x << " " << y << " " << p << endl;
                                        A[m*npixels+(outpixel2[1]-1)*mxdim+(outpixel2[0]-1)] = p;
                                    }
                                }
                            }
                        }
                    }
                }
            }
#ifdef DEBUG
            cout << "Ending interpolation" << endl;
#endif
        }
    }

    if (saveMaps) {
        vector<double> sum;
        sum.resize(npixels);

        fitsfile * mapFits;
        for (long long i=0; i<nmaps; i++) {
            double obtMin = int(intervals[0].Start());
            double obtMax = int(intervals[intervals.Count()-1].Stop());
            char out[40960] = "!";
            if(nmaps > 1) {
                if (outfile[0]=='!')
                    /*sprintf(out+1, "TA%d_TB%d_E%d-%d_TH%d-%d_SI-%g_%s",
                            int(obtMin), int(obtMax),
                            int(maps[i].emin), int(maps[i].emax),
                            int(maps[i].fovradmin), int(maps[i].fovradmax),
                            maps[i].index, outfile+1);*/
                    sprintf(out+1, "E%d-%d_TH%d-%d_SI-%g_%s",
                            int(maps[i].emin), int(maps[i].emax),
                            int(maps[i].fovradmin), int(maps[i].fovradmax),
                            maps[i].index, outfile+1);
                else
                    /*sprintf(out, "TA%d_TB%d_E%d-%d_TH%d-%d_SI-%g_%s",
                            int(obtMin), int(obtMax),
                            int(maps[i].emin), int(maps[i].emax),
                            int(maps[i].fovradmin), int(maps[i].fovradmax),
                            maps[i].index, outfile);*/
                    sprintf(out, "E%d-%d_TH%d-%d_SI-%g_%s",
                            int(maps[i].emin), int(maps[i].emax),
                            int(maps[i].fovradmin), int(maps[i].fovradmax),
                            maps[i].index, outfile);
            } else {
                if (outfile[0]=='!')
                    sprintf(out+1, "%s", outfile+1);
                else
                    sprintf(out, "%s", outfile);
            }
            //string out = maps[i].outfile;

            cout << "Creating file " << out << endl;
            if (fits_create_file(&mapFits, out, &status) != 0) {
                cerr << "Error opening file " << out << endl;
                return status;
            }
            cout << "Created file " << out << endl;

            cout << "Creating image in " << out << endl;
            int bitpix = DOUBLE_IMG;
            long naxis = 2;
            long naxes[2] = { mxdim, mxdim };
            long startpixel[2] = {1,1};
            fits_create_img(mapFits, bitpix, naxis, naxes, &status);
            // Write exp sum from all the intervals
            for (unsigned int j=0; j<npixels; j++)
                sum[j] = 0;
            for (int intvIndex=0; intvIndex<intervals.Count(); intvIndex++)
                for (unsigned int j=0; j<npixels; j++)
                    sum[j] += exposures[intvIndex][i * npixels + j];
            fits_write_pix(mapFits, TDOUBLE, startpixel, npixels, &sum[0], &status);

            double fovradmin = maps[i].fovradmin;
            double fovradmax = maps[i].fovradmax;
            double emin = maps[i].emin;
            double emax = maps[i].emax;
            double spectral_index = maps[i].index;

            fits_update_key(mapFits, TDOUBLE, "CRVAL1", &la, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "CRVAL2", &ba, NULL, &status);
            char projstr1[FLEN_FILENAME];
            char projstr2[FLEN_FILENAME];
            if (proj == AIT) {
                strcpy(projstr1,"GLON-AIT");
                strcpy(projstr2,"GLAT-AIT");
            } else {
                strcpy(projstr1,"GLON-ARC");
                strcpy(projstr2,"GLAT-ARC");
            }
            fits_update_key(mapFits, TSTRING, "CTYPE1", projstr1, NULL, &status);
            fits_update_key(mapFits, TSTRING, "CTYPE2", projstr2, NULL, &status);

            double xx = mdim/mres/2+0.5;
            fits_update_key(mapFits, TDOUBLE, "CRPIX1", &xx, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "CRPIX2", &xx, NULL, &status);
            xx = -mres;
            fits_update_key(mapFits, TDOUBLE, "CDELT1", &xx, NULL, &status);
            xx = mres;
            fits_update_key(mapFits, TDOUBLE, "CDELT2", &xx, NULL, &status);

            char unit[] = "deg";
            fits_update_key(mapFits, TSTRING, "CUNIT1", unit, NULL, &status);
            fits_update_key(mapFits, TSTRING, "CUNIT2", unit, NULL, &status);

            char strValue[72];
            strValue[68] = 0;
            strncpy(strValue, AfterSlashPtr(sarFileName), 68);
            fits_update_key(mapFits, TSTRING, "SARFILE", strValue, NULL, &status);
            strncpy(strValue, AfterSlashPtr(edpFileName), 68);
            fits_update_key(mapFits, TSTRING, "EDPFILE", strValue, NULL, &status);

            char str3[] =  "FK5";
            fits_update_key(mapFits, TSTRING, "RADESYS", str3, NULL, &status);

            xx = 2000.;
            fits_update_key(mapFits, TDOUBLE, "EQUINOX", &xx, NULL, &status);

            fits_update_key(mapFits, TDOUBLE, "LONPOLE", &lonpole, NULL, &status);
            WriteTime(mapFits, tmin, tmax);

            fits_update_key(mapFits, TDOUBLE, "MINENG", &emin, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "MAXENG", &emax, NULL, &status);
            char str1[] =  "AGILE";
            fits_update_key(mapFits, TSTRING, "TELESCOP", str1, NULL, &status);

            char str2[] =  "GRID";
            fits_update_key(mapFits, TSTRING, "INSTRUME", str2, NULL, &status);

            char str6[] =  "T";

            fits_update_key(mapFits, TSTRING, "PIXCENT", str6, NULL, &status);

            char str7[FLEN_FILENAME];
            strcpy(str7, "cm**2 s sr");
            fits_update_key(mapFits, TDOUBLE, "YTOL", &y_tol, "Pointing direction tolerance (deg)", &status);
            fits_update_key(mapFits, TDOUBLE, "ROLTOL", &roll_tol, "Roll angle tolerance (deg)", &status);
            fits_update_key(mapFits, TDOUBLE, "EARTOL", &earth_tol, "Pointing direction tolerance (deg)", &status);
            fits_update_key(mapFits, TINT, "STEP", &binstep, "Map interpolation step size", &status);
            fits_update_key(mapFits, TINT, "TIMESTEP", &timestep, "Log file step size", &status);
            fits_update_key(mapFits, TSTRING, "BUNIT", str7, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "FOV", &fovradmax, "Radius of field of view (deg)", &status);
            fits_update_key(mapFits, TDOUBLE, "FOVMIN", &fovradmin, "Minimum off-axis angle (deg)", &status);
            fits_update_key(mapFits, TDOUBLE, "ALBEDO", &albrad, "Earth zenith angle (deg)", &status);
            fits_update_key(mapFits, TINT, "PHASECOD", &phasecode, "Orbital phase code", &status);
            fits_update_key(mapFits, TDOUBLE, "INDEX", &spectral_index, NULL, &status);

            fits_update_key(mapFits, TDOUBLE, "TSTART", &obtMin, "[OBT]first event time", &status);
            fits_update_key(mapFits, TDOUBLE, "TSTOP", &obtMax, "[OBT]last event time", &status);

            char timeSys[] = "TT";
            char timeUnit[] = "s";
            double zero = 0;
            fits_update_key(mapFits, TSTRING, "TIMESYS", timeSys, "TT", &status);
            fits_update_key(mapFits, TSTRING, "TIMEUNIT", timeUnit, "TT", &status);
            fits_update_key(mapFits, TDOUBLE, "TZERO", &zero, timeSys, &status);

            fits_update_key(mapFits, TDOUBLE, "SC-Z-LII", &lp0, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "SC-Z-BII", &bp0, NULL, &status);
            fits_update_key(mapFits, TDOUBLE, "SC-LONPL", &gp0, NULL, &status);
            cout << "Closing " << out << endl;
            fits_close_file(mapFits, &status);
        }
    }
    fits_close_file(selectionFits, &status);
    fits_close_file(templateFits, &status);

    return status;
}

int EvalCounts(const char *outfile, const char *projection, double tmin,
               double tmax, double mdim, double mres, double la, double ba,
               double lonpole, double emin, double emax, double fovradmax,
               double fovradmin, double albrad, int phasecode, int filtercode,
               const char *selectionFilename,  const char *templateFilename,
               Intervals &intervals, vector< vector<int> > &counts, bool saveMaps)
{
    int status = 0;

    long mxdim = long(mdim / mres + 0.1); // dimension (in pixels) of the map
#ifdef DEBUG
    cout << "mdim: " << mdim << " mres: " << mres << " mxdim: " << mxdim << endl;
#endif
    long npixels = mxdim * mxdim;
    counts.resize(intervals.Count());
    for (int i=0; i < intervals.Count(); i++) {
        counts[i].resize(npixels);
        for (int j=0; j < npixels; j++)
            counts[i][j] = 0;
    }

    struct prjprm *prj = new(prjprm);
    prjini(prj);
    ProjectionType proj;
    if (!LoadProjection(projection, proj)) {
        cerr << "Error loading projection '" << projection << "'" << endl;
        return false;
    }

    int hdutype;
    fitsfile* selectionFits;
    if (fits_open_file(&selectionFits, selectionFilename, READONLY, &status)) {
        cerr << "ERROR opening selection file " << selectionFilename << endl;
        return -1;
    }
    fits_movabs_hdu(selectionFits, 2, &hdutype, &status);

    fitsfile* templateFits;
    if (fits_open_file(&templateFits, templateFilename, READWRITE, &status)) {
        cerr << "ERROR opening template file " << templateFilename << endl;
        return -1;
    }
    fits_movabs_hdu(templateFits, 2, &hdutype, &status);
    long oldnrows;
    fits_get_num_rows(templateFits, &oldnrows, &status);
    fits_delete_rows(templateFits, 1, oldnrows, &status);

#ifdef DEBUG
    cout << "Evaluating counts.." << endl;
#endif

    int totalCounts = 0;
    for (int intvIndex = 0; intvIndex < intervals.Count(); intvIndex++) {
#ifdef DEBUG
        cout << "Interval #" << intvIndex << endl;
#endif
        Intervals sIntv;
        sIntv.Add(intervals[intvIndex]);
        string selExpr = selection::TimesExprString(sIntv);
#ifdef DEBUG
        cout << "selExpr: " << selExpr << endl;
#endif
        fits_select_rows(selectionFits, templateFits, (char*)selExpr.c_str(), &status);
#ifdef DEBUG
        cout << "Rows from " << tempFilename << " selected" << endl;
#endif
        long nrows;
        fits_get_num_rows(templateFits, &nrows, &status);
#ifdef DEBUG
        cout << "Reading all " << nrows << " rows" << endl;
#endif

        int raColumn, decColumn;
        fits_get_colnum(templateFits, 1, (char*)"RA", &raColumn, &status);
        fits_get_colnum(templateFits, 1, (char*)"DEC", &decColumn, &status);

        double ra, dec, l, b, the, x, y, i = 0, ii = 0;
        double baa = ba * DEG2RAD;
        double laa = la * DEG2RAD;
        for (long k = 0; k<nrows; k++) {
            fits_read_col(templateFits, TDOUBLE, raColumn, k+1, 1, 1, NULL, &ra, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, decColumn, k+1, 1, 1, NULL, &dec, NULL, &status);
            Euler(ra, dec, &l, &b, 1);
            l *= DEG2RAD;
            b *= DEG2RAD;
            the = sin(b)*sin(baa)+cos(b)*cos(baa)*cos(l-laa);
            if (the < -1.0)
                the = M_PI;
            else if (the > 1.0)
                the = 0.0;
            else
                the = acos(the);

            if (proj == ARC) {
                x = RAD2DEG/Sinaa(the) * cos(b)*sin(l-laa);
                y = RAD2DEG/Sinaa(the) * (sin(b)*cos(baa) - cos(b)*sin(baa)*cos(l-laa));

                i = (int)floor(((-x+(mdim/2.))/mres));
                ii = (int)floor(((y+(mdim/2.))/mres));
            }
            else if (proj == AIT) {
                l=l-laa;

                if (l < M_PI)
                l=-l;
                else
                l=2*M_PI -l;

                x=RAD2DEG*(sqrt(2.0)*2.0*cos(b)*sin(l/2.0))/sqrt(1.0 + cos(b)*cos(l/2.0) );
                y=RAD2DEG*(sqrt(2.0)*sin(b))/sqrt(1.0 + cos(b)*cos(l/2.0) );

                i=(int)floor(((x+(mdim/2.))/mres));
                ii=(int)floor(((y+(mdim/2.))/mres));
            }

            if (inmap(i, ii, mxdim)) {
                counts[intvIndex][ii*mxdim+i]++;
                totalCounts++;
            }
        }
        if (nrows > 0)
            fits_delete_rows(templateFits, 1, nrows, &status);
    }
#ifdef DEBUG
    cout << "Ending evaluation" << endl;
#endif

    if (saveMaps) {
        vector<unsigned short> sum;
        sum.resize(npixels);
        for (unsigned int i=0; i<mxdim; i++)
            for (unsigned int j=0; j<mxdim; j++)
                sum[i*mxdim+j] = 0;
        for (int intvIndex=0; intvIndex<intervals.Count(); intvIndex++)
            for (unsigned int i=0; i<mxdim; i++)
                for (unsigned int j=0; j<mxdim; j++)
                    sum[i*mxdim+j] += counts[intvIndex][i*mxdim+j];

        int bitpix = USHORT_IMG;
        long naxis = 2;
        long naxes[2] = { mxdim, mxdim };
        fitsfile * mapFits;
        cout << "Creating file " << outfile << endl;
        if (fits_create_file(&mapFits, outfile, &status) != 0) {
            cerr << "Error opening file " << outfile << endl;
            return status;
        }
        cout << "creating Counts Map...................................." << endl;
        fits_create_img(mapFits, bitpix, naxis, naxes, &status);
        cout << "writing Counts Map with " << totalCounts << " events" << endl;
        fits_write_img(mapFits, bitpix, 1, npixels, &sum[0], &status);
        cout << "writing header........................................" << endl << endl;
        fits_update_key(mapFits, TDOUBLE, "CRVAL1", &la, NULL, &status);
        fits_update_key(mapFits, TDOUBLE, "CRVAL2", &ba, NULL, &status);
        char projstr1[FLEN_FILENAME];
        char projstr2[FLEN_FILENAME];
        if (proj == AIT) {
            strcpy(projstr1,"GLON-AIT");
            strcpy(projstr2,"GLAT-AIT");
        } else {
            strcpy(projstr1,"GLON-ARC");
            strcpy(projstr2,"GLAT-ARC");
        }
        fits_update_key(mapFits, TSTRING, "CTYPE1", projstr1, NULL, &status);
        fits_update_key(mapFits, TSTRING, "CTYPE2", projstr2, NULL, &status);
        double xx = mdim/mres/2+0.5;
        fits_update_key(mapFits, TDOUBLE, "CRPIX1", &xx, NULL, &status);
        fits_update_key(mapFits, TDOUBLE, "CRPIX2", &xx, NULL, &status);
        xx = -mres;
        fits_update_key(mapFits, TDOUBLE, "CDELT1", &xx, NULL, &status);
        xx = mres;
        fits_update_key(mapFits, TDOUBLE, "CDELT2", &xx, NULL, &status);
        char unit[] = "deg";
        fits_update_key(mapFits, TSTRING, "CUNIT1", unit, NULL, &status);
        fits_update_key(mapFits, TSTRING, "CUNIT2", unit, NULL, &status);
        char str3[] = "FK5";
        fits_update_key(mapFits, TSTRING,  "RADESYS", str3, NULL, &status);
        xx = 2000.0;
        fits_update_key(mapFits, TDOUBLE,  "EQUINOX", &xx, NULL, &status);
        fits_update_key(mapFits, TDOUBLE,  "LONPOLE", &lonpole, NULL, &status);
        WriteTime(mapFits, tmin, tmax);
        fits_update_key(mapFits, TDOUBLE,  "MINENG", &emin, NULL, &status);
        fits_update_key(mapFits, TDOUBLE,  "MAXENG", &emax, NULL, &status);
        char str1[] = "AGILE";
        fits_update_key(mapFits, TSTRING,  "TELESCOP", str1, NULL, &status);
        char str2[] = "GRID";
        fits_update_key(mapFits, TSTRING,  "INSTRUME", str2, NULL, &status);
        char str6[] = "T";
        fits_update_key(mapFits, TSTRING,  "PIXCENT", str6, NULL, &status);
        char str7[FLEN_FILENAME] = "";
        fits_update_key(mapFits, TSTRING,  "BUNIT", str7, NULL, &status);
        fits_update_key(mapFits, TDOUBLE,  "FOV", &fovradmax, "Maximum off-axis angle (deg)", &status);
        fits_update_key(mapFits, TDOUBLE,  "FOVMIN", &fovradmin, "Minimum off-axis angle (deg)", &status);
        fits_update_key(mapFits, TDOUBLE,  "ALBEDO", &albrad, "Earth zenith angle (deg)", &status);
        fits_update_key(mapFits, TINT,  "PHASECOD", &phasecode, "Orbital phase code", &status);
        fits_update_key(mapFits, TINT,  "FILTERCO", &filtercode, "Event filter code", &status);
        cout << "Closing " << outfile << endl;
        fits_close_file(mapFits, &status);
    }
    fits_close_file(selectionFits, &status);
    fits_close_file(templateFits, &status);

    return status;
}

int EvalGasMap(AgileMap &gasMap, AgileMap &expMap, const char* loresdiffuseFilename,
               const char* hiresdiffuseFilename)
{
    int status = 0;

    long outpixel[2];
    int bitpix   =  DOUBLE_IMG;
    double pixels = 0.0;
    double nulval = 1.0/0.0;
    int anynul = 0;

    double lox0 = 0, loy0 = 0;
    int lonaxis;
    long lonaxes[3];
    double lol0 = 0, lob0 = 0;
    double lodl = 0, lodb = 0;
    long lonx = 0, lony = 0;


    int hiresstatus = 0;
    double hix0 = 0, hiy0 = 0;
    int hinaxis;
    long hinaxes[3];
    double hil0 = 0, hib0 = 0;
    double hidl = 0, hidb = 0;
    long hinx = 0, hiny = 0;

    fitsfile* lodiffuseFits;
    fitsfile* hidiffuseFits;

    if (fits_open_image(&lodiffuseFits, loresdiffuseFilename, READONLY, &status) != 0) {
        printf("Errore in apertura file '%s'\n", loresdiffuseFilename);
        return status;
        }

    if (fits_open_image(&hidiffuseFits, hiresdiffuseFilename, READONLY, &hiresstatus) != 0) {
        printf("Cannot read high res diffuse file '%s'\n", hiresdiffuseFilename);
    } else {
        fits_get_img_param(hidiffuseFits, 3, &bitpix, &hinaxis, hinaxes, &hiresstatus);
        fits_read_key(hidiffuseFits,TLONG,"NAXIS1",&hinx,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TLONG,"NAXIS2",&hiny,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TDOUBLE,"CRVAL1",&hil0,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TDOUBLE,"CDELT1",&hidl,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TDOUBLE,"CRPIX1",&hix0,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TDOUBLE,"CRVAL2",&hib0,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TDOUBLE,"CDELT2",&hidb,NULL,&hiresstatus);
        fits_read_key(hidiffuseFits,TDOUBLE,"CRPIX2",&hiy0,NULL,&hiresstatus);
        if (hiresstatus != 0)
            printf("ERROR: '%s' will not be used because a FITS keyword is missing\n", hiresdiffuseFilename);
    }

    gasMap = expMap;

    fits_get_img_param(lodiffuseFits, 3, &bitpix, &lonaxis, lonaxes, &status);
    fits_read_key(lodiffuseFits,TLONG,"NAXIS1",&lonx,NULL,&status);
    fits_read_key(lodiffuseFits,TLONG,"NAXIS2",&lony,NULL,&status);
    fits_read_key(lodiffuseFits,TDOUBLE,"CRVAL1",&lol0,NULL,&status);
    fits_read_key(lodiffuseFits,TDOUBLE,"CDELT1",&lodl,NULL,&status);
    fits_read_key(lodiffuseFits,TDOUBLE,"CRPIX1",&lox0,NULL,&status);
    fits_read_key(lodiffuseFits,TDOUBLE,"CRVAL2",&lob0,NULL,&status);
    fits_read_key(lodiffuseFits,TDOUBLE,"CDELT2",&lodb,NULL,&status);
    fits_read_key(lodiffuseFits,TDOUBLE,"CRPIX2",&loy0,NULL,&status);

    long nrows = expMap.Rows();
    long ncols = expMap.Cols();

    long lodiffusepixel[3] = {1,1,1};
    long hidiffusepixel[4] = {1,1,1,1};

    bool inhiresmap = FALSE;
    double ll, bb;

    for (outpixel[0]=1;outpixel[0]<=nrows;outpixel[0]++) {
        for(outpixel[1]=1;outpixel[1]<=ncols;outpixel[1]++) {
            ll = expMap.l(outpixel[0]-1,outpixel[1]-1);
            bb = expMap.b(outpixel[0]-1,outpixel[1]-1);
            if (hiresstatus != 0)
                inhiresmap = FALSE;
            else {
                hidiffusepixel[0] = (long)(0.5 + fmod(hix0 + (ll - hil0) / hidl + 720.0/fabs(hidl), 360.0/fabs(hidl)));
                hidiffusepixel[1] = (long)(0.5 + hiy0 + (bb - hib0) / hidb);
                inhiresmap = hidiffusepixel[0] >= 1 && hidiffusepixel[0] <= hinx && hidiffusepixel[1] >= 1 && hidiffusepixel[1] <= hiny;
            }
            if (inhiresmap) {
                fits_read_pix(hidiffuseFits, TDOUBLE, hidiffusepixel, 1, &nulval, &pixels, &anynul, &hiresstatus);
                if (hiresstatus == 0 && anynul == 0) {
                    gasMap(outpixel[0]-1, outpixel[1]-1) = pixels;
                }
                else {
                    hiresstatus = 0;
                    inhiresmap = FALSE;
                }
            }
            if (!inhiresmap) {
                lodiffusepixel[0] = (long)(0.5 + fmod(lox0 + (ll - lol0) / lodl + 720.0/fabs(lodl), 360.0/fabs(lodl)));
                lodiffusepixel[1] = (long)(0.5 + loy0 + (bb - lob0) / lodb);
                if (lodiffusepixel[1] > lony)
                    lodiffusepixel[1] = lony;
                else if (lodiffusepixel[1] < 1)
                    lodiffusepixel[1] = 1;

                fits_read_pix(lodiffuseFits, TDOUBLE, lodiffusepixel, 1, NULL, &pixels, NULL, &status);
                gasMap(outpixel[0]-1, outpixel[1]-1) = pixels;
            }
            if (status) throw;
        }
    }
    fits_close_file(lodiffuseFits, &status);
    fits_close_file(hidiffuseFits, &hiresstatus);

    return status;
}

int EvalGas(const char* outfile, const char* expfile, const char* diffusefile,
            const char* hiresdiffusefile)
{
    AgileMap expMap;
    int status = expMap.Read(expfile);
    if (status) return status;

    AgileMap gasMap;
    status = EvalGasMap(gasMap, expMap, diffusefile, hiresdiffusefile);
    if (status) return status;

    std::string outstr(outfile);
    bool gz = false;
    if (outstr.size() >= 3 && outstr.compare(outstr.size() - 3, 3, ".gz") == 0) {
        gz = true;
        outstr[outstr.size()-3] = '\0';
    }
    status = gasMap.WriteWithAllMetadata(outstr.c_str());
    if (status) return status;

    FitsFile outf;
    if (!outf.Open(outstr.c_str(), READWRITE)) {
        cerr << "ERROR " << outf.Status() << " opening " << outstr << endl;
    return outf.Status();
    }
    FitsFile expf;
    if (!expf.Open(expfile)) {
        cerr << "ERROR " << expf.Status() << " opening " << expfile << endl;
    return outf.Status();
    }
    FitsFile lof;
    if (!lof.Open(diffusefile)) {
        cerr << "ERROR " << lof.Status() << " opening " << diffusefile << endl;
    return lof.Status();
    }

    std::string diff(diffusefile);
    size_t idx = diff.find_last_of("\\/");
    if (std::string::npos != idx)
        diff.erase(0, idx+1);
    std::string hidiff(hiresdiffusefile);
    idx = hidiff.find_last_of("\\/");
    if (std::string::npos != idx)
        hidiff.erase(0, idx+1);

    outf.UpdateKey("BUNIT", "(cm**2 s sr)**(-1)");
    outf.UpdateKey("EXTNAME", "GAS");
    outf.UpdateKey("SKYL", diff.c_str());
    outf.UpdateKey("SKYH", hidiff.c_str());
    char str[1024];
    if (lof.ReadKey("DH_CONF_", str, ""))
        outf.UpdateKey("DH_CONF_", str);
    if (lof.ReadKey("STAN_CON", str, ""))
        outf.UpdateKey("STAN_CON", str);
    if (lof.ReadKey("FILE_ID", str, ""))
        outf.UpdateKey("FILE_ID", str);
    char sarfile[1024];
    expf.ReadKey("SARFILE", sarfile, "");
    if (std::string(sarfile).compare("")) {
        std::string tmp(sarfile);
        std::string::size_type pos1 = sizeof("AG_GRID_")-1; // after this prefix
        std::string::size_type pos2 = tmp.find("_", pos1)+1;
        std::string::size_type pos3 = tmp.find("_", pos2)+1;
        std::string::size_type pos4 = tmp.find(".", pos3);
        if (pos2 != std::string::npos && pos3 != std::string::npos && pos4 != std::string::npos) {
            std::string dhconf = tmp.substr(pos1, pos2-pos1-1);
            std::string stancon = tmp.substr(pos2, pos3-pos2-1);
            std::string fileid = tmp.substr(pos3, pos4-pos3);
            outf.UpdateKey("DH_CONF_", dhconf.c_str());
            outf.UpdateKey("STAN_CON", stancon.c_str());
            outf.UpdateKey("FILE_ID", fileid.c_str());
        }
    }
    outf.Close();

    if (gz) {
        std::string cmd = "gzip ";
        cmd += outstr;
        system(cmd.c_str());
    }
    return status;
}

int EvalCountsInRadius(const char *outfile, double tmin,
               double tmax, double radius, double la, double ba,
               double lonpole, double emin, double emax, double fovradmax,
               double fovradmin, double albrad, int phasecode, int filtercode,
               const char *selectionFilename,  const char *templateFilename,
               Intervals &intervals, vector<int> &counts)
{
    int status = 0;

    counts.resize(intervals.Count());
    for (int i=0; i < intervals.Count(); i++) {
            counts[i] = 0;
    }

    int hdutype;
    fitsfile* selectionFits;
    if (fits_open_file(&selectionFits, selectionFilename, READONLY, &status)) {
        cerr << "ERROR opening selection file " << selectionFilename << endl;
        return -1;
    }
    fits_movabs_hdu(selectionFits, 2, &hdutype, &status);

    fitsfile* templateFits;
    if (fits_open_file(&templateFits, templateFilename, READWRITE, &status)) {
        cerr << "ERROR opening template file " << templateFilename << endl;
        return -1;
    }
    fits_movabs_hdu(templateFits, 2, &hdutype, &status);
    long oldnrows;
    fits_get_num_rows(templateFits, &oldnrows, &status);
    fits_delete_rows(templateFits, 1, oldnrows, &status);

#ifdef DEBUG
    cout << "Evaluating counts.." << endl;
#endif

    int totalCounts = 0;
    for (int intvIndex = 0; intvIndex < intervals.Count(); intvIndex++) {
#ifdef DEBUG
        cout << "Interval #" << intvIndex << endl;
#endif
        Intervals sIntv;
        sIntv.Add(intervals[intvIndex]);
        //string selExpr = selection::TimesExprString(sIntv);
        string selExpr = selection::EvtExprString(sIntv, emin, emax, albrad, fovradmax, fovradmin, phasecode, filtercode);
#ifdef DEBUG
        cout << "selExpr: " << selExpr << endl;
#endif
		cout << "selExpr: " << selExpr << endl;
        fits_select_rows(selectionFits, templateFits, (char*)selExpr.c_str(), &status);
#ifdef DEBUG
        cout << "Rows from " << tempFilename << " selected" << endl;
#endif
        long nrows;
        fits_get_num_rows(templateFits, &nrows, &status);
#ifdef DEBUG
        cout << "Reading all " << nrows << " rows" << endl;
#endif

        int raColumn, decColumn;
        fits_get_colnum(templateFits, 1, (char*)"RA", &raColumn, &status);
        fits_get_colnum(templateFits, 1, (char*)"DEC", &decColumn, &status);
        
        int enColumn, pheColumn, thetaColumn, phaseColumn, timeColumn;
        fits_get_colnum(templateFits, 1, (char*)"ENERGY", &enColumn, &status);
        fits_get_colnum(templateFits, 1, (char*)"PH_EARTH", &pheColumn, &status);
        fits_get_colnum(templateFits, 1, (char*)"THETA", &thetaColumn, &status);
        fits_get_colnum(templateFits, 1, (char*)"PHASE", &phaseColumn, &status);
		fits_get_colnum(templateFits, 1, (char*)"TIME", &timeColumn, &status);
		
        double ra, dec, l, b, the;
        double baa = ba;// * DEG2RAD;
        double laa = la;// * DEG2RAD;
        double timec, energyc, ph_earthc, thetac, phasec;
        for (long k = 0; k<nrows; k++) {
            fits_read_col(templateFits, TDOUBLE, raColumn, k+1, 1, 1, NULL, &ra, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, decColumn, k+1, 1, 1, NULL, &dec, NULL, &status);
            Euler(ra, dec, &l, &b, 1);
            //l *= DEG2RAD;
            //b *= DEG2RAD;
            fits_read_col(templateFits, TDOUBLE, enColumn, k+1, 1, 1, NULL, &energyc, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, timeColumn, k+1, 1, 1, NULL, &timec, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, pheColumn, k+1, 1, 1, NULL, &ph_earthc, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, thetaColumn, k+1, 1, 1, NULL, &thetac, NULL, &status);
            fits_read_col(templateFits, TDOUBLE, phaseColumn, k+1, 1, 1, NULL, &phasec, NULL, &status);
            
            
            double the = SphDistDeg(l, b, laa, baa);
			
			if (the < radius) {
            	totalCounts++;
            	cout << timec << " " << l << " " << b << " " << energyc << " " << thetac << " " << ph_earthc << " " << the << endl;
            }
        }
        counts[intvIndex]=totalCounts;
        if (nrows > 0)
            fits_delete_rows(templateFits, 1, nrows, &status);
    }
#ifdef DEBUG
    cout << "Ending evaluation" << endl;
#endif

    fits_close_file(selectionFits, &status);
    fits_close_file(templateFits, &status);

    return status;
}

}
