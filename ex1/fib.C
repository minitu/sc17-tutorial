#include "fib.decl.h"

#define THRESHOLD 3 // calculate sequentially if below

/* readonly */ int n; // fibonacci number
/* readonly */ CProxy_Main mainProxy;

class Main : public CBase_Main {
  double start_time;

  public:
  Main(CkArgMsg* m) {
    // determine fibonacci number
    n = 20;
    if (m->argc == 2)
      n = atoi(m->argv[1]);

    CkPrintf("\n[Exercise 1. Fibonacci Sequence]\n");

    // create chares and start computation
    start_time = CkWallTimer();
    CProxy_Fib::ckNew(n, true, CProxy_Fib());
  }

  void done() {
    CkPrintf("Elapsed time: %lf seconds\n\n", CkWallTimer() - start_time);
    CkExit();
  }
};

class Fib : public CBase_Fib {
  Fib_SDAG_CODE
  CProxy_Fib parent;
  bool is_root;

  public:
  Fib(int n, bool is_root_, CProxy_Fib parent_)
    : parent(parent_), is_root(is_root_) {
    calc(n);
  }

  int seqFib(int n) {
    return (n < 2) ? n : seqFib(n - 1) + seqFib(n - 2);
  }

  void respond(int val) {
    if (!is_root) {
      parent.response(val);
      delete this;
    }
    else {
      CkPrintf("Fibonacci number is: %d\n", val);
      mainProxy.done();
    }
  }
};

#include "fib.def.h"
