/// ///////////////////////////////////////////
///
/// HealpixMap class for writing a HEALPix FITS file
///
/// 
///
///
///
///




#ifndef _HEALPIX_MAP_
#define _HEALPIX_MAP_


enum Order { RING, NESTED };
const char* const OrderName[2] = { "RING", "NESTED" };


int Gal2Hpx(int nside, Order order, double lDeg, double bDeg);

void Hpx2Gal(int nside, Order order, int index, double* lDeg, double* bDeg);



class HealpixMap
{
public:
	HealpixMap(int nside=1, Order order=RING); /// Create an empty map
	~HealpixMap() { delete[] m_values; }

public:
	/// Table rows will contain rowOrder^2 values
	bool Save(const char* fileName, int rowOrder=1);
	bool Load(const char* fileName);

public:
	int Count() const { return 12*m_nside*m_nside; }

	int GetNSide() const { return m_nside; }
	Order GetOrder() const { return m_order; }
	
	void SetNSide(int nside);
	void SetOrder(Order order) { m_order = order; }
	
	float Val(int index) const { return m_values[index]; }
	float& Val(int index) { return m_values[index]; }
	void AddVal(int index, float value) { m_values[index]+=value; }
	
	int Gal2Hpx(double lDeg, double bDeg) const { return ::Gal2Hpx(m_nside, m_order, lDeg, bDeg); }
	void Hpx2Gal(int index, double* lDeg, double* bDeg) const { return ::Hpx2Gal(m_nside, m_order, index, lDeg, bDeg); }

private:
	int    m_nside;
	Order  m_order;
	float* m_values;
};

#endif

