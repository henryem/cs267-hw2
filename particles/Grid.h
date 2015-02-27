#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <memory>

#include "SimpleIterator.h"
#include "common.h"

class ParticleIterator;

class GridNeighborSquaresIterator;

class Grid {
private:
  friend class ParticleIterator;
  friend class GridNeighborSquaresIterator;
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
public:
  Grid(double side_length, int num_squares_per_side);
  Grid(double side_length, int num_squares_per_side, int num_particles, particle_t* particles);
  ~Grid();
  /* An iterator over particles in the same grid-square as @p or a bordering
   * square.
   */
  std::unique_ptr<SimpleIterator<particle_t&> > neighbor_iterator(const particle_t& p) const;
};

#endif /* GRID_H_ */
