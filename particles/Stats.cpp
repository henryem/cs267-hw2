/*
 * Stats.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: henrym
 */

#include <math.h>
#include "Stats.h"

void Stats::add_left(double data) {
  min = fmin(min, data);
  avg += (data - avg) / (n+1.0); // Knuth's algorithm.
  n ++;
}

void Stats::aggregate_left(Stats& addend) {
  if (addend.n > 0) {
    min = fmin(min, addend.min);
    avg += (addend.avg - avg) / (n/addend.n + 1.0);
    n += addend.n;
  }
}

Stats* Stats::aggregate(Stats& addend) {
  Stats* left_copy = new Stats(addend);
  left_copy->aggregate_left(addend);
  return left_copy;
}
