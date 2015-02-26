#include "Grid.h"

SCENARIO( "Adding and moving items in a grid", "[grid]" ) {
  GIVEN("An empty grid") {
    WHEN("we iterate over particles near some particle") {
      THEN("there are none") {
        Grid g(1.0, 1, 0, nullptr);
        std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(0.5, 0.5, 0.0, 0.0, 0.0, 0.0));
        while (i->hasNext()) {
          FAIL("there was an unexpected element in the grid");
        }
      }
    }
  }

  GIVEN("A 1x1 grid with a few elements") {
    WHEN("we iterate over particles near some particle") {
      THEN("we see all the particles") {

      }
    }
  }

  GIVEN("A 2x2 grid with a few elements, and no elements in one square") {
    WHEN("we iterate over particles near some particle") {
      THEN("we see all the particles") {

      }
    }
  }

  GIVEN("A 3x3 grid with a few elements in each square") {
    WHEN("we iterate over particles near a top-left particle") {
      THEN("we see only particles in the top 2x2 subgrid") {

      }
    }
  }
}
