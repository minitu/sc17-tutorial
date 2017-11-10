#ifndef PARTICLE_H_
#define PARTICLE_H_

// Particle object with x & y coordinates
struct Particle {
  double x, y;
  int color;

  Particle() {}
  Particle(double a, double b, int c) { x = a; y = b; color = c; }
  void pup(PUP::er& p) { p|x; p|y; p|color; }
};

#endif
