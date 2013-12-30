#ifndef _INTERVALS_
#define _INTERVALS_


class Interval
{
public:
	Interval(): m_start(0), m_stop(0) {}
	Interval(double start, double stop) { Set(start, stop); }
	~Interval() {}

	Interval& operator=(const Interval& interval) { m_start=interval.m_start; m_stop=interval.m_stop; return *this; }

	bool operator<(const Interval& intv) const { return m_stop<intv.m_start;  } /// completely preceding with no overlap
	bool operator>(const Interval& intv) const { return m_start>intv.m_stop;  } /// completely following with no overlap

	bool operator<=(const Interval& intv) const { return m_start<=intv.m_start; } /// it may overlap
	bool operator>=(const Interval& intv) const { return m_start>=intv.m_start; } /// it may overlap

	bool operator==(const Interval& intv) const { return m_start==intv.m_start && m_stop==intv.m_stop; }
	bool operator!=(const Interval& intv) const { return m_start!=intv.m_start || m_stop!=intv.m_stop; }

	friend bool Disjointed(const Interval& int1, const Interval& int2) { return int1<int2 || int2<int1; }
	friend bool Overlap(const Interval& int1, const Interval& int2) { return !Disjointed(int1, int2); }

	/// If two intervals don't overlap, Union is an interval covering both, Intersection is [0..0]
	friend Interval Union(const Interval& int1, const Interval& int2);
	friend Interval Intersection(const Interval& int1, const Interval& int2);

	void Set(double start, double stop) { if (start<stop) { m_start=start; m_stop=stop; } else { m_start=stop; m_stop=start; } }

	bool Contains(double value) const { return value>=m_start && value<=m_stop; }
	double Start() const { return m_start; }
	double Stop() const { return m_stop; }

private:
	double m_start;
	double m_stop;
};


class Intervals
{
// friend class Interval;
public:
	Intervals(): m_count(0), m_phyCount(16), m_sorted(true), m_bounds() { m_intervals = new Interval[m_phyCount]; }
	Intervals(int count): m_count(0), m_phyCount(count), m_sorted(false), m_bounds() { m_intervals = new Interval[m_phyCount]; }
	Intervals(const Intervals& other) { Copy(other); }
	~Intervals() { delete[] m_intervals; }

	Intervals& operator=(const Intervals& other) { if (this!=&other) { delete[] m_intervals; Copy(other); } return *this; }

	const Interval& operator[](int i) const { return m_intervals[i]; }
	Interval& operator[](int i) { return m_intervals[i]; }

	int Count() const { return m_count; }
	double Min() const { return m_bounds.Start(); }
	double Max() const { return m_bounds.Stop(); }

	friend Intervals Intersection(const Intervals& intervals, const Interval& interval);

	void Add(const Interval& interval);
	void Add(const Intervals& intervals);
	void Sort();
	bool Contains(double value) const;
	int IndexOf(double value) const;

private:
	int  m_count;
	int  m_phyCount;
	bool m_sorted;
	Interval m_bounds;
	Interval* m_intervals;
	void Enlarge();
	void Copy(const Intervals& other);
};


Intervals ReadIntervals(const char* intervalsFileName);


#endif

