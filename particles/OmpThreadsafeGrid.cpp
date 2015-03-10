/*
 * OmpThreadsafeGrid.cpp
 *
 *  Created on: Mar 2, 2015
 *      Author: henrym
 */

#include "omp.h"
#include "OmpThreadsafeGrid.h"

OmpThreadsafeGrid::OmpThreadsafeGrid(double side_length, int num_squares_per_side):
    Grid(side_length, num_squares_per_side),
    locks(std::vector<omp_lock_t>(num_squares_per_side*num_squares_per_side)) {
  for (unsigned int i = 0; i < locks.size(); i++) {
    omp_init_lock(&locks[i]);
  }
}

OmpThreadsafeGrid::OmpThreadsafeGrid(double side_length, int num_squares_per_side, std::vector<particle_t>& particles):
    Grid(side_length, num_squares_per_side, particles),
    locks(std::vector<omp_lock_t>(num_squares_per_side*num_squares_per_side)) {
  for (unsigned int i = 0; i < locks.size(); i++) {
    omp_init_lock(&locks[i]);
  }
}

OmpThreadsafeGrid::~OmpThreadsafeGrid() {
  for (unsigned int i = 0; i < locks.size(); i++) {
    omp_destroy_lock(&locks[i]);
  }
}

//TODO: Document thread-safety guarantee.  It is weak.
void OmpThreadsafeGrid::add(particle_t& p) {
  int square_idx = flat_idx(p);
  //FIXME: Do we need to flush before this?  When does the implied flush happen?
  //FIXME: Should technically use the scoped lock pattern here instead, in
  // case super.add() throws an exception.
  omp_set_lock(&locks[square_idx]);
  Grid::add(p);
  omp_unset_lock(&locks[square_idx]);
}

//TODO: Document thread-safety guarantee.  It is weak.
bool OmpThreadsafeGrid::remove(const particle_t& p) {
  int particle_idx = flat_idx(p);
  //FIXME: Do we need to flush before this?  When does the implied flush happen?
  //FIXME: Should technically use the scoped lock pattern here instead, in
  // case super.add() throws an exception.
  omp_set_lock(&locks[particle_idx]);
  bool result = Grid::remove(p);
  omp_unset_lock(&locks[particle_idx]);
  return result;
}
