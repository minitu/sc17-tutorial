#ifndef PARTICLE_H_
#define PARTICLE_H_

// Particle object with x & y coordinates
struct Particle {
  double x, y;

  Particle() {}
  Particle(double a, double b) { x = a; y = b; }
  void pup(PUP::er& p) { p|x; p|y; }
};

#endif
