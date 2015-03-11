/*
 * OmpThreadsafeGrid.h
 *
 *  Created on: Mar 2, 2015
 *      Author: henrym
 */

#ifndef OMPTHREADSAFEGRID_H_
#define OMPTHREADSAFEGRID_H_

#include <vector>
#include "common.h"
#include "Grid.h"

/* A Grid with a threadsafe add() operation.  The guarantee is fairly minimal;
 * it is as follows: Concurrent calls to add() are serializable.  The user is
 * responsible for flushing before any read operations and for using barriers
 * to ensure that any open iterators are destroyed before writes happen. */
class OmpThreadsafeGrid: public Grid {
private:
  /* One lock per grid square, indexed in the same way (column major).
   * The lock is held when the corresponding grid square is being modified.
   */
  std::vector<omp_lock_t> locks;
public:
  OmpThreadsafeGrid(double side_length, int num_squares_per_side);
  OmpThreadsafeGrid(double side_length, int num_squares_per_side, std::vector<particle_t>& particles);
  virtual ~OmpThreadsafeGrid();
  virtual void add(particle_t& p);
  virtual bool remove(const particle_t& p);
};

#endif /* OMPTHREADSAFEGRID_H_ */
