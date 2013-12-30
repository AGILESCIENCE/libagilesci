


#ifndef _PIL_PARAMS_
#define _PIL_PARAMS_




enum PilType { PilNone, PilBool, PilInt, PilReal, PilString };

struct PilDescription
{
	PilType     paramType;
	const char* paramName;
	const char* paramPrompt;
};


class PilValue
{
public:
	PilValue(): m_type(PilNone) {}
	PilValue(int value): m_type(PilNone) { Set(value); }
	PilValue(double value): m_type(PilNone) { Set(value); }
	PilValue(float value): m_type(PilNone) { Set(value); }
	PilValue(const char* value): m_type(PilNone) { Set(value); }
	PilValue(const PilValue& pilValue): m_type(PilNone) { Copy(pilValue); }
	~PilValue() { Clear(); }

	PilValue& operator=(const PilValue& pilValue) { if (this!=&pilValue) Copy(pilValue); return *this; }

	/// Value info
	PilType GetType() const { return m_type; }

	operator bool() const { return boolVal; }
	operator int() const { return intVal; }
	operator double() const { return floatVal; }
	operator const char*() const { return stringValue; }

	bool Get(bool& value) const { value=boolVal; return m_type==PilBool; }
	bool Get(int& value) const { value=intVal; return m_type==PilInt; }
	bool Get(double& value) const { value=floatVal; return m_type==PilReal; }
	bool Get(const char*& value) const { value=stringValue; return m_type==PilString; }

	int GetBool() const { return boolVal; }
	int GetInt() const { return intVal; }
	double GetReal() const { return floatVal; }
	const char* GetStr() const { return stringValue; }

	/// Setting the value
	void SetType(PilType pilType); /// Make it of type PilNone, or set the default value
	void Set(bool value) { Clear(); m_type=PilBool; boolVal=value; }
	void Set(int value) { Clear(); m_type=PilInt; intVal=value; }
	void Set(double value) { Clear(); m_type=PilReal; floatVal=value; }
	void Set(const char* value);

	PilValue& operator=(bool value) { Set(value); return *this; }
	PilValue& operator=(int value) { Set(value); return *this; }
	PilValue& operator=(double value) { Set(value); return *this; }
	PilValue& operator=(const char* value) { Set(value); return *this; }

private:
	PilType m_type;
	union {
		bool   boolVal;
		int    intVal;
		double floatVal;
		char*  stringValue;
		};
	void Clear() { if (m_type==PilString) delete[] stringValue; m_type=PilNone; }
	void SetStr(const char* str);
	void Copy(const PilValue& pilValue);
};

/*
void MakeNewFormat(const char* const paramNames[], const PilType* paramTypes, const char* const paramDesc[]);
*/

class PilParams
{
public:
	PilParams(const PilDescription params[]);
	~PilParams() { for (int i=0; i<=m_count; ++i) delete[] m_names[i]; delete[] m_names; }

	/// Load function to handle the errors yourself.
	/// If an error occurs accessing a parameter its name is given by parName
	bool Load(int argC, char* argV[], int* pilErr, const char*& parName);
	/// Load function that prints the errors to standard error
	bool Load(int argC, char* argV[]) { const char* s; return Load(argC, argV, 0, s); }
	/// Both Load functions return true on success

	/// Printing the parameters and their values (after loading)
	/// The original descriptor is assumed if params==0
	/// If printEmpty is true, the parameters that were not loaded are printed too
	void Print(const PilDescription params[]=0, bool printEmpty = false) const;

	int ParfileCount() const { return m_parfileCount; } /// Number of params on file
	int Count() const { return m_count; } /// Number of params by construction
	PilType GetType(const char* name) const { return m_values[Index(name)].GetType(); }

	/// Getting or setting parameters values
	const PilValue& GetValue(const char* name) const { return m_values[Index(name)]; }
	PilValue& GetValue(const char* name) { return m_values[Index(name)]; }

	const PilValue& operator[](const char* name) const { return m_values[Index(name)]; }
	PilValue& operator[](const char* name) { return m_values[Index(name)]; }

	bool GetValue(const char* name, bool& value) const { const PilValue& val = GetValue(name); return val.Get(value); }
	bool GetValue(const char* name, int& value) const { const PilValue& val = GetValue(name); return val.Get(value); }
	bool GetValue(const char* name, double& value) const { const PilValue& val = GetValue(name); return val.Get(value); }
	bool GetValue(const char* name, char* value) const;

	bool GetBoolValue(const char* name) const { const PilValue& val = GetValue(name); return val.GetBool(); }
	int GetIntValue(const char* name) const { const PilValue& val = GetValue(name); return val.GetInt(); }
	double GetRealValue(const char* name) const { const PilValue& val = GetValue(name); return val.GetReal(); }
	const char* GetStrValue(const char* name) const { const PilValue& val = GetValue(name); return val.GetStr(); }

private:
	const PilDescription* m_params; /// Stores the parameters given at construction
	int       m_count; /// Number of parameters given at construction
	char**    m_names;
	PilValue* m_values;
	int       m_parfileCount; /// Number of parameters actually in the command line
	int Index(const char* name) const;
	void Print(int index, bool printEmpty) const;
};





#endif
