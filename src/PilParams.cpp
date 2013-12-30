


#include <cstring>
#include <iostream>

#include "pil.h"
#include "PilParams.h"
#include "MathUtils.h"



using namespace std;


static char* StrCopy(const char* str)
{
int len = strlen(str);
char* stringCopy = new char[len+1];
strcpy(stringCopy, str);
return stringCopy;
}




void PilValue::SetType(PilType pilType)
{
Clear();
switch (pilType) {
	case PilBool: Set(false); break;
	case PilInt: Set(0); break;
	case PilReal: Set(0.0); break;
	case PilString: Set("");
	}
}


void PilValue::Set(const char* str)
{
Clear();
m_type = PilString;
stringValue = StrCopy(str);
}


void PilValue::Copy(const PilValue& pilValue)
{
Clear();
switch (pilValue.m_type) {
	case PilBool: Set(pilValue.boolVal); break;
	case PilInt: Set(pilValue.intVal); break;
	case PilReal: Set(pilValue.floatVal); break;
	case PilString: Set(pilValue.stringValue);
	}
}



/**
PilParams::PilParams(const char* const paramNames[], const PilType* paramTypes)
{
for (m_count=0; paramNames[m_count]; ++m_count) ;
m_names = new char*[m_count+1];
m_values = new PilValue[m_count+1];
m_names[m_count] = 0;
m_values[m_count] = PilValue();	/// For missing values

for (int i=0; i<m_count; ++i) {
	m_names[i] = StrCopy(paramNames[i]);
	m_values[i].SetType(paramTypes[i]);
	}
}
*/


/** Seems unused
static PilType PilTypeOf(const char* strType)
{
switch (strType[0]) {
	case 'B': return PilBool;
	case 'I': return PilInt;
	case 'R': return PilReal;
	case 'S': return PilString;
	}
return PilNone;
}
*/


PilParams::PilParams(const PilDescription params[])
{
m_params = params;
for (m_count=0; params[m_count].paramType!=PilNone; ++m_count) ;
m_names = new char*[m_count+1];
m_values = new PilValue[m_count+1];
m_names[m_count] = 0;
m_values[m_count] = PilValue();	/// For missing values

for (int i=0; i<m_count; ++i) {
	m_names[i] = StrCopy(params[i].paramName);
	m_values[i].SetType(params[i].paramType);
	}
m_parfileCount = 0;
}

int PilParams::Index(const char* name) const
{
for (int i=0; i<m_count; ++i)
	if (strcmp(name, m_names[i])==0)
		return i;
return m_count;	/// For missing names
}



bool PilParams::Load(int argC, char* argV[], int* pilErr, const char*& parName)
{
if (pilErr) {
	parName = "";
	*pilErr = PIL_OK;
	}
int status = PILInit(argC, argV);
if (status!=PIL_OK) {
	if (pilErr)
		*pilErr = status;
	else
		cerr << "PILInit error " << status << ": " << PIL_err_table[PIL_ERR_BASE-status] << endl;
		// cerr << "PILInit error " << status << endl;
	return false;
	}

status = PILGetNumParameters(&m_parfileCount);
if (status!=PIL_OK) {
	if (pilErr)
		*pilErr = status;
	else
		cerr << "PILGetNumParameters error " << status << ": " << PIL_err_table[PIL_ERR_BASE-status] << endl;
	return false;
	}

int count = m_count;
if (m_parfileCount<m_count)
	count = m_parfileCount;

int intVal;
double realVal;
char strVal[4096];
for (int i=0; i<count; ++i) {
	switch (m_values[i].GetType()) {
	case PilNone:
		m_values[i].SetType(PilNone);
		break;
	case PilBool:
		status = PILGetBool(m_names[i], &intVal);
		if (status==PIL_OK)
			m_values[i].Set(bool(intVal));
		break;
	case PilInt:
		status = PILGetInt(m_names[i], &intVal);
		if (status==PIL_OK)
			m_values[i].Set(intVal);
		break;
	case PilReal:
		status = PILGetReal(m_names[i], &realVal);
		if (status==PIL_OK)
			m_values[i].Set(realVal);
		break;
	case PilString:
		status = PILGetString(m_names[i], strVal);
		if (status==PIL_OK)
			m_values[i].Set(strVal);
		}
	if (status!=PIL_OK) {
		if (pilErr) {
			*pilErr = status;
			parName = m_names[i];
			}
		else
			cerr << "PIL Error " << status << " (" << PIL_err_table[PIL_ERR_BASE-status] << ") loading parameter " << m_names[i] << endl;
		return false;
		}
	}
return true;
}


/**
static void PrintNiceDouble(double val)
{
char dblStr[256];
sprintf(dblStr, "%f", val);
char* point = strchr(dblStr, '.');
if (strchr(dblStr, '.')) {
	int len = strlen(dblStr);
	while (len>1 && dblStr[len-1]=='0')
		--len;
	if (len>1 && dblStr[len-1]=='.')
		--len;
	dblStr[len] = 0;
	}
cout << dblStr << endl;
}
*/

void PilParams::Print(int index, bool printEmpty) const
{
switch (m_values[index].GetType()) {
case PilNone:
	if (printEmpty)
		cout << "[Empty]" << endl;
	break;
case PilBool:
	cout << (m_values[index].GetBool() ? "Yes" : "No") << endl;
	break;
case PilInt:
	cout << m_values[index].GetInt() << endl;
	break;
case PilReal:
	cout << DoubleStr(m_values[index]) << endl;
	/// PrintNiceDouble(m_values[index]);
	break;
case PilString:
	cout << m_values[index].GetStr() << endl;
	}
}


/*
void PilParams::Print(const char* const paramDesc[], bool printEmpty) const
{
for (int i=0; i<m_count; ++i) {
	cout << paramDesc[i] << ": ";
	Print(i, printEmpty);
	}
}
*/

void PilParams::Print(const PilDescription params[], bool printEmpty) const
{
if (!params)
	params = m_params;
for (int i=0; i<m_count; ++i) {
	if (strcmp(params[i].paramName, m_names[i]))
		cout << "[ERROR: " << m_names[i] << " expected] ";
	cout << params[i].paramPrompt << ": ";
	Print(i, printEmpty);
	}
}

/**
void MakeNewFormat(const char* const paramNames[], const PilType* paramTypes, const char* const paramDesc[])
{
cout << "const PilDescription c_params[] = {" << endl;
int i=0;
while (paramTypes[i]) {
	cout << "\t{ ";
	switch (paramTypes[i]) {
		case PilBool: cout << "PilBool"; break;
		case PilInt: cout << "PilInt"; break;
		case PilReal: cout << "PilReal"; break;
		case PilString: cout << "PilString";
		}
 	cout << ", \"" <<paramNames[i] << "\", \"" << paramDesc[i] << "\" }," << endl;
	++i;
	}

cout << "\t{ PilNone, \"\", \"\" }};" << endl;
}
*/