// Exercise4: Navier-Stokes equations
//
// Goal: solve ther Navier-Stokes equations
//
// Reference: see https://goo.gl/JK8EUI for a detailed description of possible 
//            expression terms.
//
// Compile-and-run: 
// > cd build
// > make exercise4
// > cd ..
// > build/exercise4 init/exercise4.dat.2d
//

#include "AMDiS.h"

using namespace AMDiS;

// inflow velocity profile:
struct G : AbstractFunction<double, WorldVector<double> >
{
  G(int comp, double vel = 1.0) 
    : comp(comp), vel(vel) 
  {}
  
  double operator()(WorldVector<double> const& x) const 
  {
    // return the velocity profile
  }
private:
  int comp;   // the space component
  double vel; // the maximum inflow velocity
};

// solve: dt(u) - laplace(u) = f(x) in Omega,    u = 0 on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("ns");
  prob.initialize(INIT_ALL);
  
  ProblemInstat instat("ns", prob);
  instat.initialize(INIT_ALL);

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt", prob.getNumComponents());
  
  int dow = Global::getGeo(WORLD);
  double* tau = adaptInfo.getTimestepPtr();
  
  // shortcuts for solution vector
  WorldVector<DOFVector<double>*> u;
  for (int i = 0; i < dow; ++i) {
    u[i] = prob.getSolution(i);
  }
  
  for (int i = 0; i < dow; ++i) {
    // add operators here
  }
  
  for (int i = 0; i < dow; ++i) {
    // add boundary conditions here
  }
  
  // ===== set initial solution =====
  for (int i = 0; i < dow; ++i)
    // set initial solution component
  
  AdaptInstationary adaptInstat("adapt", prob, adaptInfo, instat, adaptInfo);
  adaptInstat.adapt();

  AMDiS::finalize();
}
