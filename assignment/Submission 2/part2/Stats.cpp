/*
 * Stats.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: henrym
 */

#include <math.h>
#include "Stats.h"

Stats::~Stats() { }

void Stats::add_left(double data) {
  min = fmin(min, data);
  avg += (data - avg) / (n+1.0); // Knuth's algorithm.
  n ++;
}

void Stats::aggregate_left(const Stats& addend) {
  if (addend.n > 0) {
    min = fmin(min, addend.min);
    avg += (addend.avg - avg) / (n/addend.n + 1.0);
    n += addend.n;
  }
}

const Stats& Stats::aggregate(const Stats& addend) const {
  Stats* left_copy = new Stats(addend);
  left_copy->aggregate_left(addend);
  return *left_copy;
}

const Stats& Stats::operator+(const Stats& addend) {
  return aggregate(addend);
}
