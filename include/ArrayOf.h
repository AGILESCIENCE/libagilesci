




#ifndef _ARRAY_OF_
#define _ARRAY_OF_

// Date: 2012-07-21

#include <stdarg.h>


template<int N, class T> class ArrayOf
{
public:
	ArrayOf(): m_size(0), m_values(0) { for (int i=0; i<=N; ++i) m_dims[i] = 0; }
	ArrayOf(int dim1, ...);
	ArrayOf(const int dimList[]);
	ArrayOf(const ArrayOf& m) { Copy(m); }
	~ArrayOf() { delete[] m_values; }

	/// Resizing
	void ResizeTo(int dim1, ...);
	void ResizeTo(const int dimList[]);
	void ResizeTo(const ArrayOf& m) { ResizeTo(m.m_dims); }

	/// General information
	int Size() const { return m_size; }
	int Dims() const { return N; }
	int Dim(int dim) const { return m_dims[dim]; }
	int DimOf(int index, int dim) const; 		/// Get the index along a given dimension
	void DimOf(int index, int* dimArr) const;	/// Get all indices

	/// For two-dimensional arrays
	int Rows() const { return m_dims[0]; } /// Alias for two dimensions matrices
	int Cols() const { return m_dims[1]; } /// Alias for two dimensions matrices
	int RowOf(int index) const { return index/m_dims[1]; }
	int ColOf(int index) const { return index%m_dims[1]; }
	void Transpose();
	void TransposeTo(ArrayOf& dest) const; /// dest cannot be self (use Transpose() in that case)
	/// See also ArrayOf operator*(const ArrayOf& a, const ArrayOf& b);

	T First() const { return m_size ? m_values[0] : T(0); }
	T Last() const { return m_size ? m_values[m_size-1] : T(0); }

	/// Find indices in arrays of increasing values
	bool IsSorted() const;
	ArrayOf& Sort();
	int LeftIndex(T value) const;	/// Return index in [-1..Size()-1] so that value>=arr[index]
	int RightIndex(T value) const;	/// Return index in [0..Size()] so that value<=arr[index]
	int LinearIndex(T value) const;	/// Closest index if the values increase linearly
	int GeomIndex(T value) const;	/// Closest index if the values increase geometrically

	/// Global assignment
	ArrayOf& operator=(const ArrayOf& m) { if (this!=&m && N==m.Dims()) { delete[] m_values; Copy(m); } return *this; }
	ArrayOf& operator=(T value) { return SetElements(value); }
	ArrayOf& operator=(T* values) { return SetElements(values); }

	ArrayOf& SetElements(T value) { for (int i=0; i<m_size; ++i) m_values[i] = value; return *this; }
	ArrayOf& SetElements(T* values) { for (int i=0; i<m_size; ++i) m_values[i] = values[i]; return *this; }

	typedef bool (*ConstArrayFunction)(int index, const T& value);
	bool Iterate(ConstArrayFunction arrFunc) const;

	/// Operations with variable number of indices

	int AbsIndex(int index1, int index2) const { return index1*m_dims[1]+index2; }
	int AbsIndex(int index1, int index2, int index3, ...) const;

	T operator()(int index) const { return m_values[index]; }
	T& operator()(int index) { return m_values[index]; }

	T operator()(int index1, int index2) const { return m_values[index1*m_dims[1]+index2]; }
	T& operator()(int index1, int index2) { return m_values[index1*m_dims[1]+index2]; }

	T operator()(int index1, int index2, int index3, ...) const;
	T& operator()(int index1, int index2, int index3, ...);

	/// Direct access
	operator const T*() const { return m_values; }
	operator T*() { return m_values; }
	const T* Buffer() const { return m_values; }
	T* Buffer() { return m_values; }
	T operator[](int linearIndex) const { return m_values[linearIndex]; }
	T& operator[](int linearIndex) { return m_values[linearIndex]; }

	/// Math operators
	ArrayOf& operator+=(const ArrayOf& m) { for (int i=0; i<m_size; ++i) m_values[i]+=m.m_values[i]; return *this; }
	ArrayOf& operator-=(const ArrayOf& m) { for (int i=0; i<m_size; ++i) m_values[i]-=m.m_values[i]; return *this; }

	ArrayOf& operator*=(T scalar) { for (int i=0; i<m_size; ++i) m_values[i] *= scalar; return *this; }
	ArrayOf& operator/=(T scalar) { for (int i=0; i<m_size; ++i) m_values[i] /= scalar; return *this; }

	friend ArrayOf operator*(const ArrayOf& m, T scalar) { ArrayOf m2(m); m2*=scalar; return m2; }
	friend ArrayOf operator/(const ArrayOf& m, T scalar) { ArrayOf m2(m); m2/=scalar; return m2; }

private:
	int  m_dims[N+1];	/// Size of each dimension
	int  m_size;		/// Number of elements (product of all m_dims)
	T*   m_values;	
	void Copy(const ArrayOf& m);
	void CopyDims(const int dimList[]);
	void Delete() { delete[] m_values; }
	int OutIndex(const T& value) const;
	void MergeSort(const ArrayOf& arr1, const ArrayOf& arr2);
};


/// Basic vectors
typedef ArrayOf<1, int>    VecI;
typedef ArrayOf<1, float>  VecF;
typedef ArrayOf<1, double> VecD;

/// Basic matrices
typedef ArrayOf<2, int>    MatI;
typedef ArrayOf<2, float>  MatF;
typedef ArrayOf<2, double> MatD;

/// Basic 3D matrices
typedef ArrayOf<3, int>    Mat3I;
typedef ArrayOf<3, float>  Mat3F;
typedef ArrayOf<3, double> Mat3D;



template<int N, class T> ArrayOf<N, T> operator*(const ArrayOf<N, T>& a, const ArrayOf<N, T>& b)
{
if (a.Dims()!=2 || b.Dims()!=2 || a.Dim(1)!=b.Dim(0))
	return ArrayOf<2, T>();
ArrayOf<2, T> result(a.Dim(0), b.Dim(1));
result = T(0);
for (int i=0; i<a.Dim(0); ++i)
	for (int j=0; j<b.Dim(1); ++j)
		for (int k=0; k<a.Dim(1); ++k)
			result(i,j) += a(i,k)*b(k,j);
return result;
}


template<int N, class T> bool ArrayOf<N, T>::Iterate(ConstArrayFunction arrFunc) const
{
bool ok = true;
for (int i=0; i<m_size && ok; ++i)
	ok = arrFunc(i, m_values[i]);
return ok;
}



template<int N, class T> void ArrayOf<N, T>::Copy(const ArrayOf<N, T>& m)
{
for (int i=0; i<=N; ++i)
	m_dims[i] = m.m_dims[i];
m_size = m.m_size;
m_values = new T[m_size];
for (int j=0; j<m_size; ++j)
	m_values[j] = m.m_values[j];
}



template<int N, class T> ArrayOf<N, T>::ArrayOf(int dim1, ...)
{
m_dims[N] = 0;
m_dims[0] = dim1;
m_size = dim1;
va_list vl;
va_start(vl, dim1);
for (int i=1; i<N; ++i) {
	int dimN = va_arg(vl, int);
	m_dims[i] = dimN;
	m_size *= dimN;
	}
va_end(vl);
m_values = new T[m_size];
}




template<int N, class T> ArrayOf<N, T>::ArrayOf(const int dimList[])
{
m_size = 1;
for (int i=0; i<N; ++i)
	m_size *= m_dims[i] = dimList[i];
m_dims[N] = 0;
m_values = new T[m_size];
}


template<int N, class T> int ArrayOf<N, T>::DimOf(int index, int dim) const
{
if (dim>=N)
	return 0;
int mult = m_size;
for (int i=0; i<dim; ++i) {
	mult /= m_dims[i];
	index %= mult;
	}
mult /= m_dims[dim];
return index/mult;
}


template<int N, class T> void ArrayOf<N, T>::DimOf(int index, int* dimArr) const
{
int mult = m_size;
for (int i=0; i<N; ++i) {
	mult /= m_dims[i];
	dimArr[i] = index/mult;
	index %= mult;
	}
}


template<int N, class T> void ArrayOf<N, T>::ResizeTo(const int dimList[])
{
int size = m_size;
m_size = 1;
for (int i=0; i<N; ++i)
	m_size *= m_dims[i] = dimList[i];
if (size!=m_size) {
	T* values = new T[m_size];
	int minSize = m_size<size?m_size:size;
	for (int i=0; i<minSize; ++i)
		values[i] = m_values[i];
	delete[] m_values;
	m_values = values;
	for (; size<m_size; ++size)
		m_values[size] = T(0);
	}
}


template<int N, class T> void ArrayOf<N, T>::ResizeTo(int dim1, ...)
{
int dims[N];
dims[0] = dim1;
if (N>1) {
	va_list vl;
	va_start(vl, dim1);
	for (int i=1; i<N; ++i)
		dims[i] = va_arg(vl, int);
	va_end(vl);
	}
ResizeTo(dims);
}


template<int N, class T> void ArrayOf<N, T>::TransposeTo(ArrayOf& dest) const
{
if (this==&dest || N!=2 || dest.Dims()!=2)
	return;
dest.ResizeTo(Cols(), Rows());
const ArrayOf& mat(*this);
for (int i=0; i<Rows(); ++i)
	for (int j=0; j<Cols(); ++j)
		dest(j,i) = mat(i,j);
}

template<int N, class T> void ArrayOf<N, T>::Transpose()
{
if (Cols()==Rows()) {
	ArrayOf& mat(*this);
	T swap;
	for (int i=0; i<Rows()-1; ++i)
		for (int j=i+1; j<Cols(); ++j) {
			swap = mat(i, j);
			mat(i, j) = mat(j, i);
			mat(j, i) = swap;
			}
	}
else {
	ArrayOf mat;
	TransposeTo(mat);
	for (int i=0; i<=N; ++i)
		m_dims[i] = mat.m_dims[i];
	delete m_values;
	m_values = mat.m_values;
	mat.m_values = 0;
	}
}




#ifdef EvalArrayOfAbsIndex
#error EvalArrayOfAbsIndex should not be defined at this point
#endif

#define EvalArrayOfAbsIndex \
	int mult = m_size/m_dims[0]; \
	int index = index1*mult; \
	mult /= m_dims[1]; \
	index += index2*mult; \
	mult /= m_dims[2]; \
	index += index3*mult; \
	if (N>3) { \
		va_list vl; \
		va_start(vl, index3); \
		for (int i=3; i<N; ++i) { \
			mult /= m_dims[i]; \
			int indexN = va_arg(vl, int); \
			index += indexN*mult; \
			} \
		va_end(vl); \
		}

/*
#define ArrayOfBinIndices \
	int mult = m_size/m_dims[0]; \
	int index = index1*mult; \
	mult /= m_dims[1]; \
	index += index2*mult; \
	mult /= m_dims[2]; \
	index += index3*mult; \
	if (N>3) { \
		va_list vl; \
		va_start(vl, index3); \
		for (int i=3; i<N; ++i) { \
			mult /= m_dims[i]; \
			int indexN = va_arg(vl, int); \
			index += indexN*mult; \
			} \
		va_end(vl); \
		} \
	return m_values[index]
*/

template<int N, class T> int ArrayOf<N, T>::AbsIndex(int index1, int index2, int index3, ...) const
{
EvalArrayOfAbsIndex;
return index;
}

template<int N, class T> T ArrayOf<N, T>::operator()(int index1, int index2, int index3, ...) const
{
EvalArrayOfAbsIndex;
return m_values[index];
}

template<int N, class T> T& ArrayOf<N, T>::operator()(int index1, int index2, int index3, ...)
{
EvalArrayOfAbsIndex;
return m_values[index];
}

#undef EvalArrayOfAbsIndex




template<int N, class T> bool ArrayOf<N, T>::IsSorted() const
{
for (int i=1; i<m_size; ++i)
	if (m_values[i]<m_values[i-1])
		return false;
return true;
}


template<int N, class T> ArrayOf<N, T>& ArrayOf<N, T>::Sort()
{
if (m_size<=1)
	return *this;
else if (m_size==2) {
	if (m_values[0]>m_values[1]) {
		T temp(m_values[0]);
		m_values[0] = m_values[1];
		m_values[1] = temp;
		}
	return *this;
	}
int mid = m_size/2;
ArrayOf arr1(mid);
ArrayOf arr2(m_size-mid);
arr1.SetElements(m_values);
arr2.SetElements(m_values+mid);
arr1.Sort();
arr2.Sort();
MergeSort(arr1, arr2);
return *this;
}

template<int N, class T> void ArrayOf<N, T>::MergeSort(const ArrayOf& arr1, const ArrayOf& arr2)
{
ResizeTo(arr1.Size()+arr2.Size());
int idx1 = 0;
int idx2 = 0;
int index = 0;
while (index<m_size)
	if (idx2>=arr2.Size() || (idx1<arr1.Size() && arr1[idx1]<arr2[idx2]))
		m_values[index++] = arr1[idx1++];
	else
		m_values[index++] = arr2[idx2++];
}



template<int N, class T> int ArrayOf<N, T>::OutIndex(const T& value) const
{
if (m_size<1 || value<m_values[0])
	return -1;
else if (m_size<2)
	return 0;
else if (value>=m_values[m_size-1])
	return m_size-1;
return -2;
}


template<int N, class T> int ArrayOf<N, T>::LeftIndex(const T value) const
{
int first = OutIndex(value);
if (first>-2)
	return first;
first = 0;
int last = m_size-1;
int mid = last/2;
while (true) {
	if (value<m_values[mid])
		if (mid==first+1)
			return first;
		else
			last = mid;
	else
		if (last==mid+1)
			return mid;
		else
			first = mid;
	mid = (first+last)/2;
	}
return 0;
}


template<int N, class T> int ArrayOf<N, T>::RightIndex(T value) const
{
int index = LeftIndex(value);
if (index<=0)
	return 0;
else if (value==m_values[index])
	--index;
return index;
}


template<int N, class T> int ArrayOf<N, T>::LinearIndex(T value) const
{
int index = LeftIndex(value);
if (index<0)
	return 0;
else if (index==m_size-1)
	return m_size-1;
return (value-m_values[index]<=m_values[index+1]-value) ? index : index+1;
}


template<int N, class T> int ArrayOf<N, T>::GeomIndex(T value) const
{
int index = LeftIndex(value);
if (index<0)
	return 0;
else if (index==m_size-1)
	return m_size-1;
return (value/m_values[index]<=m_values[index+1]/value) ? index : index+1;
}


#endif
