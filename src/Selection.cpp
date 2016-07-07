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
#include <fstream>
#include <sstream>
#include <cstdio>
#include <fitsio.h>
#include <string.h>

#include <pil.h>
#include <FitsUtils.h>

#include "Selection.h"

//#define DEBUG 1

using std::string;
using std::stringstream;
using std::cout;
using std::cerr;
using std::endl;

namespace selection
{

int pil_chridx(char *p, char c)
{
    int  l;

    for (l = 0; p[l]; l++)
        if (c == p[l])
            return(l);
    return(-1);
}

int pil_curly_append(char **p, char *s, int nchars)
{
    int l, nl;
    char *np;

    if (NULL == p)  return(PIL_NUL_PTR);
    if (nchars < 0) return(PIL_BAD_ARG);
    if (0 == nchars) {
        if (NULL != *p)
            return(PIL_OK);
    }
    else {
        if (NULL == s)
            return(PIL_NUL_PTR);
    }

    if (NULL == *p)
        l = 0;
    else
        l = strlen(*p);

    nl = nchars + l + 1;      // + 1 - for EOS

    if (NULL == *p)
        np = (char *)PIL_malloc(nl);
    else
        np = (char *)PIL_realloc((void *)*p, nl);

    if (NULL == np) {
        if (NULL != *p) {
            PIL_free(*p);    // on error deallocate string
            *p = NULL;
        }
        return(PIL_NO_MEM);
    }

    if (nchars > 0)
        memcpy(np + l, s, nchars); // append new string
    np[nl - 1] = 0;          // signal EOS

    *p = np;
    return(PIL_OK);
}

int pil_curly_expand(char *src, char **dst)
{
    char ev[PIL_CURLYSIZE];
    int sl, l, l2, l3, tl, r;
    char *p;

    if (NULL == src)
        return(PIL_NUL_PTR);
    if (NULL == dst)
        return(PIL_NUL_PTR);
    *dst = NULL;

    sl = tl = 0;
    for (;;) {
        if (-1 == (l = pil_chridx(src + sl + tl, '$')))
            break;          // no '$' from current pos, so terminate
        if ('{' != src[sl + tl + l + 1]) {
            tl += l + 1;
            continue;       // '$' not followed by '{'
        }

        l += tl;            // total number of bytes _BEFORE_ '${...}'
        if (-1 == (l2 = pil_chridx(src + sl + l + 2, '}')))
            break;          // no closing '}'
        if (PIL_OK != (r = pil_curly_append(dst, src + sl, l)))
            return(r);      // copy part _BEFORE_ '$'
        sl += l + 2;        // position on 1st char after '${'

        if (l2 > 0) {       // if non-empty env.var name
            memcpy(ev, src + sl, l2);    // then copy its name to the temporary buffer
            ev[l2] = 0;     // signal end of string
            p = getenv(ev); // get environment variable
            l3 = 0;
            if (NULL != p)
                l3 = strlen(p); // ... and compute length of that variable (0 - if notdef)
            if (PIL_OK != (r = pil_curly_append(dst, p, l3)))
                return(r);  // copy value of environment variable
        }
        sl += l2 + 1;       // env.var.name + '}'
        tl = 0;
    }

    return(pil_curly_append(dst, src + sl, strlen(src + sl))); // copy remaining part
}

string TimesExprString(const Intervals& intvs)
{
    if (intvs.Count() <= 0)
        return string("");

    stringstream str(ios_base::out);
    str.precision(6);
    if (intvs.Count() == 1)
        str << "TIME >= " << std::fixed << intvs[0].Start() << " && TIME < " << intvs[0].Stop();
    else {
        str << "(";
        const char* sep = "";
        for (int i=0; i<intvs.Count(); i++) {
            str << sep << "(TIME >= " << std::fixed << intvs[i].Start() << " && TIME < " << intvs[i].Stop() << ")";
            sep = " || ";
        }
        str << ")";
    }
    return str.str();
}

string LogExprString(const Intervals& intvs, int phasecode, int timeStep)
{
    if (intvs.Count() <= 0)
        return string("");
    stringstream str(std::ios_base::out);
    str << TimesExprString(intvs);
    str << " && LIVETIME > 0 && LOG_STATUS == 0 && MODE == 2";
    if ((phasecode & 1) == 1)
        str << " && PHASE .NE. 0";
    if ((phasecode & 2) == 2)
        str << " && PHASE .NE. 1";
    if ((phasecode & 4) == 4)
        str << " && PHASE .NE. 2";
    if ((phasecode & 8) == 8)
        str << " && PHASE .NE. 3";
    if ((phasecode & 16) == 16)
        str << " && PHASE .NE. 4";
    str << " && ((#ROW == 1) || (#ROW == (#ROW/" << timeStep << ") *" << timeStep << "))";
    return str.str();
}

void CopyFile(const char* iname, const char* oname)
{
    std::ifstream ifs(iname, std::ios::binary);
    std::ofstream ofs(oname, std::ios::binary);
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
}

int MakeSelection(const char *fileList, Intervals& selection,
                  string expr, const char *selectionFilename,
                  const char *templateFilename)
{
    FILE *fp = fopen(fileList, "r");
    if (!fp) {
        cerr << "Error opening file " << fileList << endl;
        return -100;
    }

    cout << "MakeSelection intervals: " << endl << selection.String() << endl;

    char scratchFilename[FLEN_FILENAME];
    tmpnam(scratchFilename);

    FitsFile selectionFits;
    int skippedFiles = 0;
    int status = 0;
    char buffer[40960];
    bool selectionOpened = false;
    while (fgets(buffer, 40960, fp) && !status) {
        char name[FLEN_FILENAME];
        double t1=0, t2=0;
        sscanf(buffer, "%s %lf %lf", name, &t1, &t2);

        Interval logIntv(t1, t2);
        Intervals sel = Intersection(selection, logIntv);

        char* nname = 0;
        pil_curly_expand(name, &nname);
        if (sel.Count()) {
            if (skippedFiles) {
                cout << "Skipped " << skippedFiles << " files." << endl;
                skippedFiles = 0;
            }
            cout << "Selecting from: " << nname << endl;
        }
        else
            skippedFiles++;

        for (int i=0; i<sel.Count() && !status; i++) {
            Intervals selInt;
            selInt.Add(sel[i]);

            string sin(nname);
            if(sin.compare(sin.size()-3, sin.size(), ".gz") == 0) {
                string command = string("gunzip -c ") + nname + string(" > ") + scratchFilename;
#ifdef DEBUG
                cout << command << endl;
#endif
                system(command.c_str());
            }
            else
                CopyFile(nname, scratchFilename);
            FitsFile scratchFits;
            if (!scratchFits.Open(scratchFilename)) {
                cerr << "Failed opening scratch: " << scratchFilename << endl;
                PIL_free(nname);
                fclose(fp);
                return -100;
            }
            scratchFits.MoveAbsHDU(2);

            if (!selectionOpened) {
                if(!selectionFits.Open(selectionFilename, READWRITE)) {
                    CopyFile(scratchFilename, templateFilename);
                    FitsFile templateFits(templateFilename, READWRITE);
                    if (!templateFits.Ok()) {
                        char strerr[100];
                        templateFits.Error(strerr);
                        cerr << "Failed opening template: " << templateFilename  << " " << strerr << endl;
                        PIL_free(nname);
                        fclose(fp);
                        return -100;
                    }
                    templateFits.MoveAbsHDU(2);
                    long nrows;
                    fits_get_num_rows(templateFits, &nrows, &status);
                    fits_delete_rows(templateFits, 1, nrows, &status);
                    templateFits.Close();

                    selectionFits.Detach();
                    CopyFile(templateFilename, selectionFilename);
                    if (!selectionFits.Open(selectionFilename, READWRITE)) {
                        char strerr[100];
                        selectionFits.Error(strerr);
                        cerr << "Failed opening selection: " << selectionFilename  << " " << strerr << endl;
                        PIL_free(nname);
                        fclose(fp);
                        return -100;
                    }
                }
                selectionFits.MoveAbsHDU(2);
                selectionOpened = true;
            }
#ifdef DEBUG
            long nrows, nrows2;
            fits_get_num_rows(selectionFits, &nrows, &status);
            cout << "Log expr: " << expr << endl;
#endif
            fits_select_rows(scratchFits, selectionFits, (char*)expr.c_str(), &status);
#ifdef DEBUG
            fits_get_num_rows(selectionFits, &nrows2, &status);
            cout << nrows2-nrows << " rows total" << endl;
            cout << "Selection file total: " << nrows2 << endl;
#endif
            scratchFits.Delete();
        }
        PIL_free(nname);
    }
    fclose(fp);
    if (skippedFiles)
        cout << "Skipped " << skippedFiles << " files." << endl;

    long allnrows;
    fits_get_num_rows(selectionFits, &allnrows, &status);
    selectionFits.Close();
    cout << "Selected " << allnrows  << " rows." << endl;
    if (allnrows == 0)
        return -118;

    return status;
}

}
