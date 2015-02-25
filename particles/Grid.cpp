/*
 * Grid.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: henrym
 */

#include "Grid.h"
#include "SimpleIterator.h"
#include "common.h"



/* Accomplishes the same thing as the following Scala one-liner:
 *   grid_squares_iterator.flatmap(g.squares)
 * ...but with a bunch of complicated and specialized C++ code.
 */
class ParticleIterator : public SimpleIterator<particle_t&> {
private:
  const Grid& g;
  SimpleIterator<int>& grid_squares_iterator;
  // Uninitialized until has_started is true.  Nonsense if has_finished is
  // true.
  // Otherwise this always points to a valid square.  (If there are no
  // valid squares, hasNext() will never return true and has_started will
  // never become true.  The contract guarantees that next() will never be
  // called in that case.)
  int square_idx;
  // Uninitialized until has_started is true.  Nonsense if has_finished is
  // true.
  // The index of the next particle that will be returned.
  int particle_in_square_idx;
  bool has_started;
  bool has_finished;

public:
  ParticleIterator(const Grid& g_v, SimpleIterator<int>& grid_squares_iterator_v):
      g(g_v),
      grid_squares_iterator(grid_squares_iterator_v),
      square_idx(0),
      particle_in_square_idx(0),
      has_started(false),
      has_finished(false) { }

  ~ParticleIterator() {
    delete &grid_squares_iterator;
  }

  bool hasNext() const {
    if (has_finished) {
      return false;
    }
    if (!has_started) {
      return grid_squares_iterator.hasNext();
    }
    return true;
  }

  particle_t& next() {
    if (!has_started) {
      square_idx = grid_squares_iterator.next();
      has_started = true;
    }
    particle_t& p = *(g.squares[square_idx][particle_in_square_idx]);
    particle_in_square_idx++;
    if (particle_in_square_idx >= g.squares[square_idx].size()) {
      if (grid_squares_iterator.hasNext()) {
        square_idx = grid_squares_iterator.next();
        particle_in_square_idx = 0;
      } else {
        has_finished = true;
      }
    }
    return p;
  }
};


/* Iterates over the raw indices of the 9 squares neighboring a given grid
 * square.  Border squares are skipped.  The iteration order is arbitrary.
 * -1 is returned when the iterator is exhausted.
 */
class GridNeighborSquaresIterator : public SimpleIterator<int> {
private:
  // The following invariant is maintained: (current_square_x, current_square_y)
  // is within the bounding box defined by the following edge coordinates, OR
  // hasNext() is false.
  const int top_edge_x;
  const int left_edge_y;
  const int bottom_edge_x;
  const int right_edge_y;
  const Grid& g;
  int current_square_x;
  int current_square_y;

public:
  GridNeighborSquaresIterator(const Grid& g_v, int square_x, int square_y):
      top_edge_x(max(square_x-1, 0)),
      left_edge_y(max(square_y-1, 0)),
      bottom_edge_x(min(square_x+1, g.num_squares_per_side)),
      right_edge_y(min(square_y+1, g.num_squares_per_side)),
      g(g_v),
      current_square_x(top_edge_x),
      current_square_y(left_edge_y) { }
  ~GridNeighborSquaresIterator() { }
  int next() {
    const int current_square_idx = g.square_idx_to_flat_idx(current_square_x, current_square_y);
    current_square_x++;
    if (current_square_x > bottom_edge_x) {
      current_square_x = top_edge_x;
      current_square_y++;
    }
    return current_square_idx;
  }
  bool hasNext() const {
    return current_square_y < right_edge_y;
  }
};


Grid::Grid(double side_length_v, int num_squares_per_side_v):
    side_length(side_length_v),
    num_squares_per_side(num_squares_per_side_v),
    square_length(side_length_v / num_squares_per_side_v),
    num_squares(num_squares_per_side*num_squares_per_side),
    squares(new std::vector<particle_t*>[num_squares_per_side*num_squares_per_side]) {
}

Grid::Grid(double side_length_v, int num_squares_per_side_v, int num_particles, particle_t* particles):
    Grid::Grid(side_length_v, num_squares_per_side_v) {
  for (int particle_idx = 0; particle_idx < num_particles; particle_idx++) {
    particle_t* p = &particles[particle_idx];
    squares[flat_idx(*p)].push_back(p);
  }
}

Grid::~Grid() {
  //FIXME: Still not really sure how this works in C++.
  delete squares;
}

SimpleIterator<particle_t&>& Grid::neighbor_iterator(const particle_t& p) const {
  //FIXME: Not sure whether it's idiomatic to use the heap here.  There is a
  // memory leak.
  SimpleIterator<int>* neighbor_squares = new GridNeighborSquaresIterator(*this, square_x(p.x), square_y(p.y));
  SimpleIterator<particle_t&>* neighbor_particles = new ParticleIterator(*this, *neighbor_squares);
  return *neighbor_particles;
}
