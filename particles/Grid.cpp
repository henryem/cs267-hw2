/*
 * Grid.cpp
 *
 *  Created on: Feb 20, 2015
 *      Author: henrym
 */

#include <memory>

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
  std::unique_ptr<SimpleIterator<int> > grid_squares_iterator;
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
  unsigned int particle_in_square_idx;
  bool has_started;
  bool has_finished;

public:
  ParticleIterator(const Grid& g_v, std::unique_ptr<SimpleIterator<int> > grid_squares_iterator_v):
      g(g_v),
      grid_squares_iterator(std::move(grid_squares_iterator_v)),
      square_idx(0),
      particle_in_square_idx(0),
      has_started(false),
      has_finished(false) { }

  ~ParticleIterator() { }

  bool hasNext() {
    // has_finished is set by next().
    if (has_finished) {
      return false;
    }
    if (!has_started) {
      // Search for a non-empty grid square and start there.  Since hasNext()
      // must be called before next(), we don't have to worry about repeating
      // this initialization in next().
      while (grid_squares_iterator->hasNext()) {
        int next_square_idx = grid_squares_iterator->next();
        if (!g.square(next_square_idx).empty()) {
          this->square_idx = next_square_idx;
          has_started = true;
          return true;
        }
      }
      // No grid squares had elements.
      has_finished = true;
      return false;
    } else {
      return true;
    }
  }

  particle_t& next() {
    particle_t* p = g.square(square_idx)[particle_in_square_idx];
    particle_in_square_idx++;
    while (particle_in_square_idx >= g.square(square_idx).size()) {
      if (grid_squares_iterator->hasNext()) {
        square_idx = grid_squares_iterator->next();
        particle_in_square_idx = 0;
      } else {
        has_finished = true;
        break;
      }
    }
    return *p;
  }
};


/* Iterates over the raw indices of the 9 squares neighboring a given grid
 * square.  Border squares are skipped.  The iteration order is arbitrary.
 * -1 is returned when the iterator is exhausted.
 */
class GridNeighborSquaresIterator : public SimpleIterator<int> {
private:
  const Grid& g;
  // The following invariant is maintained: (current_square_x, current_square_y)
  // is within the bounding box defined by the following edge coordinates, OR
  // hasNext() is false.  The bounding box is top- and left-inclusive and bottom-
  // and right-exclusive.
  const int top_edge_x;
  const int left_edge_y;
  const int bottom_edge_x;
  const int right_edge_y;
  int current_square_x;
  int current_square_y;

public:
  GridNeighborSquaresIterator(const Grid& g_v, int square_x, int square_y):
      g(g_v),
      top_edge_x(max(square_x-1, 0)),
      left_edge_y(max(square_y-1, 0)),
      bottom_edge_x(min(square_x+2, g.num_squares_per_side)),
      right_edge_y(min(square_y+2, g.num_squares_per_side)),
      current_square_x(top_edge_x),
      current_square_y(left_edge_y) { }
  ~GridNeighborSquaresIterator() { }
  int next() {
    const int current_square_idx = g.square_idx_to_flat_idx(current_square_x, current_square_y);
    current_square_x++;
    if (current_square_x >= bottom_edge_x) {
      current_square_x = top_edge_x;
      current_square_y++;
    }
    return current_square_idx;
  }
  bool hasNext() {
    return current_square_y < right_edge_y;
  }
};

class RangeIterator : public SimpleIterator<int> {
private:
  int current;
  int end;
public:
  /* Inclusive on the left, exclusive on the right. */
  RangeIterator(int start, int end_v): current(start), end(end_v) { }
  int next() {
    return current++;
  }
  bool hasNext() {
    return current < end;
  }
};


Grid::Grid(double side_length_v, int num_squares_per_side_v):
    side_length(side_length_v),
    num_squares_per_side(num_squares_per_side_v),
    square_length(side_length_v / num_squares_per_side_v),
    num_squares(num_squares_per_side*num_squares_per_side),
    squares(std::unique_ptr<std::vector<Square> >(new std::vector<Square>(num_squares_per_side*num_squares_per_side))) {
}

Grid::Grid(double side_length_v, int num_squares_per_side_v, std::vector<particle_t>& particles):
    Grid::Grid(side_length_v, num_squares_per_side_v) {
  for (std::vector<particle_t>::iterator i = particles.begin(); i != particles.end(); i++) {
    particle_t& p = *i;
    add(p);
  }
}

Grid::~Grid() { }

std::unique_ptr<SimpleIterator<particle_t&> > Grid::neighbor_iterator(const particle_t& p) const {
  int x = square_x(p.x);
  int y = square_y(p.y);
  return std::unique_ptr<SimpleIterator<particle_t&> >(new ParticleIterator(
      *this,
      std::unique_ptr<SimpleIterator<int> >(new GridNeighborSquaresIterator(*this, x, y))));
}

std::unique_ptr<SimpleIterator<particle_t&> > Grid::subgrid(int subgrid_idx, int num_subgrids) const {
  //TODO: For now just returns a range of subblocks according to column-major
  // order.  This makes the surface area bigger, resulting in more
  // communication when this is used as a blocking structure.
  int num_squares_per_subgrid = div_round_up(num_squares, num_subgrids);
  int start = subgrid_idx * num_squares_per_subgrid;
  int end = min(start + num_squares_per_subgrid, num_squares);
  return std::unique_ptr<SimpleIterator<particle_t&> >(new ParticleIterator(
      *this,
      std::unique_ptr<SimpleIterator<int> >(new RangeIterator(start, end))));
}

void Grid::add(particle_t& p) {
  square(flat_idx(p)).push_back(&p);
}

bool Grid::remove(const particle_t& p) {
  Square& s = square(p);
  for (Square::iterator i = s.begin(); i != s.end(); i++) {
    particle_t& possible_p = **i;
    if (&possible_p == &p) {
      s.erase(i);
      return true;
    }
  }
  return false;
}
