#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <memory>

#include "SimpleIterator.h"
#include "common.h"

class ParticleIterator;

class GridNeighborSquaresIterator;

class RangeIterator;

class Grid {
protected:
  friend class ParticleIterator;
  friend class GridNeighborSquaresIterator;
  friend class RangeIterator;
  //TODO: Should be a vector of weak_ptr<particle_t>, but that seems a little
  // complicated, so I'm punting for now.  This class does not own the
  // particles.
  typedef std::vector<particle_t*> Square;
  const double side_length;
  const int num_squares_per_side;
  const double square_length;
  const int num_squares;
  std::unique_ptr<std::vector<Square> > squares;
  // The flat index of the grid square containing particle p.
  inline int flat_idx(const particle_t& p) const {
    return flat_idx(p.x, p.y);
  }
  inline int square_x(double x) const {
    return (int) (x / square_length);
  }
  inline int square_x(const particle_t& p) const {
    return (int) (p.x / square_length);
  }
  inline int square_y(double y) const {
    return (int) (y / square_length);
  }
  inline int square_y(const particle_t& p) const {
    return (int) (p.y / square_length);
  }
  // The flat index of the grid square containing point (x,y).  x is the distance
  // from the top edge of the grid and y is the distance from the left.
  inline int flat_idx(double x, double y) const {
    return square_idx_to_flat_idx(square_x(x), square_y(y));
  }
  inline int square_idx_to_flat_idx(int square_x, int square_y) const {
    return square_x + square_y*num_squares_per_side;
  }
  inline Square& square(int square_idx) const {
    return (*squares)[square_idx];
  }
  inline Square& square(const particle_t& p) const {
    return square(flat_idx(p));
  }
public:
  Grid(double side_length, int num_squares_per_side);
  Grid(double side_length, int num_squares_per_side, std::vector<particle_t>& particles);
  virtual ~Grid();
  /* An iterator over particles in the same grid-square as @p or a bordering
   * square.
   */
  std::unique_ptr<SimpleIterator<particle_t&> > neighbor_iterator(const particle_t& p) const;
  /* An iterator over particles in subgrid @subgrid_idx out of @num_subgrids
   * total subgrids.  It is guaranteed that the disjoint union of subgrid(0, k),
   * subgrid(1, k), ..., subgrid(k-1, k) equals the set of particles in the
   * grid.  Beyond that, subclasses may divide themselves arbitrarily but
   * should attempt to create contiguous subgrids with small average surface
   * area (e.g. O(sqrt(subgrid size))).
   */
  std::unique_ptr<SimpleIterator<particle_t&> > subgrid(int subgrid_idx, int num_subgrids) const;
  /* Add @p to the grid.  Subclasses may define different thread-safety
   * semantics for this method. */
  virtual void add(particle_t& p);
  /* Remove @p from the grid.  Returns true if @p was present in the grid, and
   * false otherwise.  Equality is defined by identity (i.e. &particle == &p).
   * Only the first instance of @p is removed, if there are multiples (though
   * it would probably be wrong to store multiple copies of a particle in a
   * grid).
   *
   * Note that this method may take O(n) time, where n is the number of
   * particles currently in the grid.
   *
   * Subclasses may define different thread-safety semantics for this method.
   */
  virtual bool remove(const particle_t& p);
};

#endif /* GRID_H_ */
