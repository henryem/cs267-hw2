#include "Grid.h"

SCENARIO( "Adding and moving items in a grid", "[grid]" ) {
  GIVEN("An empty grid") {
    Grid g(1.0, 1, std::vector<particle_t>(0));
    WHEN("we iterate over particles near some particle") {
      std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(0.5, 0.5));
      THEN("there are none") {
        REQUIRE_FALSE(i->hasNext());
      }
    }
  }

  GIVEN("A 1x1 grid with a few elements") {
    std::vector<particle_t> ps = {particle_t(0.2, 0.2), particle_t(0.3, 0.3)};
    Grid g(1.0, 1, ps);
    WHEN("we iterate over particles near some particle") {
      std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(0.5, 0.5));
      THEN("we see all the particles") {
        int num_particles_observed = 0;
        while (i->hasNext()) {
          particle_t& p = i->next();
          num_particles_observed++;
        }
        REQUIRE(num_particles_observed == 2);
      }
    }
    WHEN("a particle is modified") {
      std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(0.5, 0.5));
      ps[0].x = 0.4;
      THEN("the change is reflected in particles seen by neighbor_iterator") {
        while (i->hasNext()) {
          particle_t& p = i->next();
          REQUIRE((p.x == 0.3 || (p.x == 0.4 && p.y == 0.2)));
        }
      }
    }
  }

  GIVEN("A 2x2 grid with a few elements, and no elements in one square") {
    std::vector<particle_t> ps = {particle_t(0.2, 0.2), particle_t(0.3, 0.3), particle_t(0.8, 0.8), particle_t(0.8, 0.3)};
    Grid g(1.0, 2, ps);
    WHEN("we iterate over particles near some particle") {
      std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(0.5, 0.5));
      THEN("we see all the particles") {
        int num_particles_observed = 0;
        while (i->hasNext()) {
          particle_t& p = i->next();
          num_particles_observed++;
        }
        REQUIRE(num_particles_observed == 4);
      }
    }
  }

  GIVEN("A 3x3 grid with a few elements in each square") {
    std::vector<particle_t> ps = {
        particle_t(0.5, 0.5),
        particle_t(1.5, 0.5),
        particle_t(2.5, 0.5),
        particle_t(0.5, 1.5),
        particle_t(1.5, 1.5),
        particle_t(2.5, 1.5),
        particle_t(0.5, 2.5),
        particle_t(1.5, 2.5),
        particle_t(2.5, 2.5),
    };
    Grid g(3.0, 3, ps);
    WHEN("we iterate over particles near a top-left particle") {
      std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(0.5, 0.5));
      THEN("we see 4 particles") {
        int num_particles_observed = 0;
        while (i->hasNext()) {
          particle_t& p = i->next();
          num_particles_observed++;
        }
        REQUIRE(num_particles_observed == 4);
      }
      THEN("we see only particles in the top 2x2 subgrid") {
        while (i->hasNext()) {
          particle_t& p = i->next();
          REQUIRE((&p == &ps[0] || &p == &ps[1] || &p == &ps[3] || &p == &ps[4]));
        }
      }
    }

    GIVEN("A 5x5 grid with an element in each square") {
      std::vector<particle_t> ps = {
          particle_t(0.5, 0.5),
          particle_t(1.5, 0.5),
          particle_t(2.5, 0.5),
          particle_t(3.5, 0.5),
          particle_t(4.5, 0.5),
          particle_t(0.5, 1.5),
          particle_t(1.5, 1.5),
          particle_t(2.5, 1.5),
          particle_t(3.5, 1.5),
          particle_t(4.5, 1.5),
          particle_t(0.5, 2.5),
          particle_t(1.5, 2.5),
          particle_t(2.5, 2.5),
          particle_t(3.5, 2.5),
          particle_t(4.5, 2.5),
          particle_t(0.5, 3.5),
          particle_t(1.5, 3.5),
          particle_t(2.5, 3.5),
          particle_t(3.5, 3.5),
          particle_t(4.5, 3.5),
          particle_t(0.5, 4.5),
          particle_t(1.5, 4.5),
          particle_t(2.5, 4.5),
          particle_t(3.5, 4.5),
          particle_t(4.5, 4.5),
      };
      Grid g(5.0, 5, ps);
      WHEN("we iterate over particles near the middle") {
        std::unique_ptr<SimpleIterator<particle_t&> > i = g.neighbor_iterator(particle_t(2.5, 2.5));
        THEN("we see 9 particles") {
          int num_particles_observed = 0;
          while (i->hasNext()) {
            particle_t& p = i->next();
            num_particles_observed++;
          }
          REQUIRE(num_particles_observed == 9);
        }
      }
    }
  }
}
