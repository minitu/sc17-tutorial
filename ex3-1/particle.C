#include <stdlib.h>
#include <vector>
#include "pup_stl.h"
#include "Particle.h"
#include "particle.decl.h"

#define RANGE 1.0
#define MAX_ITERATIONS 100

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellProxy;
/* readonly */ int particles_per_cell;
/* readonly */ int cell_dimension;

class Main: public CBase_Main {
  int iteration;
  double start_time;

  public:
    Main(CkArgMsg* m) {
      if (m->argc < 3) {
        CkAbort("USAGE: ./charmrun +p<number_of_processors>"
                "./particle <number of particles per cell> <size of array>");
      }

      mainProxy = thisProxy;
      iteration = 1;
      particles_per_cell = atoi(m->argv[1]);
      cell_dimension = atoi(m->argv[2]);
      delete m;

      CkPrintf("\n[Exercise 3-1. Particles: Load Imbalance & Projections]\n");

      // create the grid and start the simulation
      start_time = CkWallTimer();
      cellProxy = CProxy_Cell::ckNew(cell_dimension, cell_dimension);
      cellProxy.run();
    }

    void printTotal(int result) {
      CkPrintf("Total: %d\n", result);
    }

    void printMax(int result) {
      CkPrintf("Maximum: %d\n", result);
    }

    void done() {
      if (iteration++ < MAX_ITERATIONS) {
        cellProxy.run();
      }
      else {
        CkPrintf("Elapsed time: %lf\n\n", CkWallTimer() - start_time);
        CkExit();
      }
    }
};

// This class represents a cell of the simulation, and each cell contains
// a vector of particles. At each iteration, the cell perturbs the particles
// and moves them to neighboring cells as necessary.
class Cell: public CBase_Cell {
  public:
    double x_min, x_max, y_min, y_max;
    int received;

    std::vector<Particle> particles;

    // particles to be transferred to neighborring cells
    std::vector<Particle> N, E, O, W, S;

    Cell() {
      initializeBounds();

      // introduce artificial load imbalance
      if (CkMyPe() == 0 || CkMyPe() == (CkNumPes() - 1)) {
        populateCell(0);
      }
      else {
        populateCell(particles_per_cell);
      }

      received = 0;
    }

    Cell(CkMigrateMessage* m) {}

    void pup(PUP::er &p) {
      CBase_Cell::pup(p);
      p|x_min; p|y_min; p|x_max; p|y_max;
      p|received;
      p|particles;
      p|N; p|E; p|O; p|W; p|S;
    }

    void run() {
      // clear previously collected particles
      N.clear(); E.clear(); O.clear(); W.clear(); S.clear();

      // move particles
      for (int i = 0; i < particles.size(); i++) {
        perturb(&particles[i]);
      }

      // update particles remaining in current cell
      particles = O;

      // send particle vectors to neighbors
      int x = thisIndex.x;
      int y = thisIndex.y;
      thisProxy(wrap(x), wrap(y + 1)).receive(N);
      thisProxy(wrap(x + 1), wrap(y)).receive(W);
      thisProxy(wrap(x - 1), wrap(y)).receive(E);
      thisProxy(wrap(x), wrap(y - 1)).receive(S);
    }

    void receive(std::vector<Particle> incoming) {
      particles.insert(particles.end(), incoming.begin(), incoming.end());

      if (++received == 4) {
        received = 0;

        // reductions for finding total and max number of particles
        reduce();

        // back to Main
        contribute(CkCallback(CkReductionTarget(Main, done), mainProxy));
      }
    }

  private:
    void reduce() {
      int num_particles = particles.size();
      contribute(sizeof(int), &num_particles, CkReduction::sum_int, CkCallback(CkReductionTarget(Main, printTotal), mainProxy));
      contribute(sizeof(int), &num_particles, CkReduction::max_int, CkCallback(CkReductionTarget(Main, printMax), mainProxy));
    }

    void initializeBounds() {
      x_min = RANGE * thisIndex.x / cell_dimension;
      x_max = RANGE * (thisIndex.x + 1) / cell_dimension;
      y_min = RANGE * thisIndex.y / cell_dimension;
      y_max = RANGE * (thisIndex.y + 1 ) / cell_dimension;
    }

    double randomWithin(double min, double max) {
      double random = drand48();
      return min + random * (max - min);
    }

    void populateCell(int initialElements) {
      for (int element = 0; element < initialElements; element++) {
        Particle p = Particle(randomWithin(x_min, x_max), randomWithin(y_min, y_max));
        particles.push_back(p);
      }
    }

    void perturb(Particle* particle) {
      float max_delta = 0.01; // determines particle speed
      double delta_x = drand48()*max_delta - max_delta/2.0;
      double delta_y = drand48()*max_delta - max_delta/2.0;

      // particle moves into either x or into y direction randomly
      double direction = drand48();
      if (direction >= 0.5)
        particle->x += delta_x;
      else
        particle->y += delta_y;

      double x = particle->x;
      double y = particle->y;
      double new_x = x;
      double new_y = y;

      // wrap around if particles go out of range
      if (new_x > RANGE) new_x -= 1.0;
      else if (new_x < 0.0) new_x += 1.0;
      if (new_y > RANGE) new_y -= 1.0;
      else if (new_y < 0.0) new_y += 1.0;

      // create a particle and put it into the correct vector to be sended
      Particle temp(new_x, new_y);
      if (y > y_max) {
        N.push_back(temp);
      } else if (x < x_min) {
        E.push_back(temp);
      } else if (x > x_max) {
        W.push_back(temp);
      } else if (y < y_min) {
        S.push_back(temp);
      } else {
        O.push_back(temp);
      }
    }

    int wrap(int w) {
      if (w >= cell_dimension) return 0;
      else if (w < 0) return cell_dimension - 1;
      else return w;
    }
};

#include "particle.def.h"
