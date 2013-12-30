


#include "Intervals.h"
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


Interval Union(const Interval& int1, const Interval& int2)
{
Interval merged(int2);
if (int1.m_start<int2.m_start)
	merged.m_start = int1.m_start;
if (int1.m_stop>int2.m_stop)
	merged.m_stop = int1.m_stop;
return merged;
}



Interval Intersection(const Interval& int1, const Interval& int2)
{
if (Disjointed(int1, int2))
	return Interval(0, 0);
Interval intersection;
intersection.m_start = int1.m_start>int2.m_start?int1.m_start:int2.m_start;
intersection.m_stop = int1.m_stop<int2.m_stop?int1.m_stop:int2.m_stop;
return intersection;
}




Intervals Intersection(const Intervals& intervals, const Interval& interval)
{
Intervals intersection;
for (int i=0; i<intervals.Count(); ++i)
	if (Overlap(intervals[i], interval))
		intersection.Add(Intersection(intervals[i], interval));
intersection.Sort();
return intersection;
}



void Intervals::Add(const Interval& interval)
{
if (m_count>=m_phyCount)
	Enlarge();
m_intervals[m_count++] = interval;
if (m_count<=1)
	m_bounds = interval;
else {
	m_sorted = false;
/**
	m_sorted = m_sorted && m_bounds.Stop()<interval.Start();
cout << (m_sorted?"Adding sorted":"Adding UNsorted") << endl;
*/
	m_bounds = Union(m_bounds, interval);
	}
}


void Intervals::Add(const Intervals& intervals)
{
int count = intervals.Count();
for (int i=0; i<count; ++i)
	Add(intervals[i]);
Sort();
}


void Intervals::Sort()
{
m_sorted = true;
if (m_count<2)
	return;

/// Sorting
for (int i=0; i<m_count-1; ++i)
	for (int j=i+1; j<m_count; ++j)
		if (m_intervals[j]<=m_intervals[i]) {
			Interval temp(m_intervals[j]);
			m_intervals[j] = m_intervals[i];
			m_intervals[i] = temp;
			}

/// Packing
Interval* intervals = new Interval[m_phyCount];
int count = 0;
for (int i=0; i<m_count; ++i) {
	if (i==m_count-1)
		intervals[count] = m_intervals[i];
	else if (m_intervals[i]<m_intervals[i+1])
		intervals[count] = m_intervals[i];
	else {
		intervals[count] = Union(m_intervals[i], m_intervals[i+1]);
		++i;
		}
	++count;
	}
delete[] m_intervals;
m_intervals = intervals;
m_count = count;
}


bool Intervals::Contains(double value) const
{
if (m_bounds.Contains(value))
	for (int i=0; i<m_count; ++i)
		if (value<=m_intervals[i].Stop())
			return value>=m_intervals[i].Start();
return false;
}

int Intervals::IndexOf(double value) const
{
if (m_bounds.Contains(value))
	for (int i=0; i<m_count; ++i)
		if (value<=m_intervals[i].Stop())
			if (value>=m_intervals[i].Start())
				return i;
return -1;
}


void Intervals::Enlarge()
{
m_phyCount += 16;
Interval* intervals = new Interval[m_phyCount];
for (int i=0; i<m_count; ++i)
	intervals[i] = m_intervals[i];
delete[] m_intervals;
m_intervals = intervals;
}

void Intervals::Copy(const Intervals& other)
{
m_phyCount = other.m_phyCount;
m_count = other.m_count;
m_sorted = other.m_sorted;
m_bounds = other.m_bounds;
m_intervals = new Interval[m_phyCount];
for (int i=0; i<m_count; ++i)
	m_intervals[i] = other.m_intervals[i];
}



static string Trim(string s)
{
int len = s.length();
if (!len)
	return s;
int a = 0;
while (a<len && s.at(a)<=32)
	++a;
if (a==len)
	return s="";
int z = len-1;
while (z>0 && s.at(z)<=32)
	--z;
return s.substr(a, z-a+1);
}

Intervals ReadIntervals(const char* intervalsFileName)
{
Intervals intervals;
int lineNum = 0;
string line;
ifstream ifile(intervalsFileName);
if (ifile.is_open()) {
	do {
		getline(ifile, line);
		line = Trim(line);
		++lineNum;
		if (line.length()!=0 && line.at(0)!='!') {
			stringstream str(line, ios_base::in);
			string name;
			double start;
			double stop;
			str >> start >> stop;
			if (stop<=start)
				cerr << "WARNING in line " << lineNum << ": start greater than stop" << endl;
			intervals.Add(Interval(start, stop));
			}
		} while (!ifile.eof());
	intervals.Sort();
	}
else
	cerr << "Error opening file " << intervalsFileName << endl;
return intervals;
}
