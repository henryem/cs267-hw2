/*
 * Stats.h
 *
 *  Created on: Feb 20, 2015
 *      Author: henrym
 */

#ifndef STATS_H_
#define STATS_H_

#include <math.h>

class Stats {
public:
  double min;
  double avg;
  int n;

  Stats(): min(INFINITY), avg(0.0), n(0) { }
	Stats(double minv, double avgv, int nv): min(minv), avg(avgv), n(nv) { }
	Stats(const Stats& original): min(original.min), avg(original.avg), n(original.n) { }
	~Stats();

	void add_left(double data);
	void aggregate_left(const Stats& addend);
	const Stats& aggregate(const Stats& addend) const;
	const Stats& operator+(const Stats& addend);
};

#endif /* STATS_H_ */
