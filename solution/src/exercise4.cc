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

struct G : AbstractFunction<double, WorldVector<double> >
{
  G(int comp, double vel = 1.0) 
    : comp(comp), vel(vel) 
  {}
  
  double operator()(WorldVector<double> const& x) const 
  {
    return comp == 0 ? 4.0*vel*x[1]*(1.0 - x[1]) : 0.0;
  }
private:
  int comp;
  double vel;
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
  double nu = 1.0, vel = 1.0;
  Parameters::get("parameters->nu", nu);
  Parameters::get("parameters->vel", vel);
  
  // shortcuts for solution vector
  WorldVector<DOFVector<double>*> u;
  for (int i = 0; i < dow; ++i) {
    u[i] = prob.getSolution(i);
  }
  
  for (int i = 0; i < dow; ++i) {
    // 1/tau * u + (u * nabla) u - nu*laplace(u)
    Operator* opL = new Operator(prob.getFeSpace(i), prob.getFeSpace(i));
    addZOT(opL, 1.0/var(tau)); // 1/tau * u
    addSOT(opL, nu);
    for (int j = 0; j < dow; ++j)
      addFOT(opL, valueOf(u[j]), j, GRD_PHI);
    prob.addMatrixOperator(opL, i, i);
    
    // 1/tau * u_old
    Operator* opRhs = new Operator(prob.getFeSpace(i));
    addZOT(opRhs, valueOf(u[i])/var(tau)); // 1/tau * u
    prob.addVectorOperator(opRhs, i);
  
    // pressure operators
    Operator* opP = new Operator(prob.getFeSpace(i), prob.getFeSpace(dow));
    addFOT(opP, 1.0, i, GRD_PSI);
    prob.addMatrixOperator(opP, i, dow); 
  
    // divergence operators
    Operator* opDiv = new Operator(prob.getFeSpace(dow), prob.getFeSpace(i));
    addFOT(opDiv, 1.0, i, GRD_PHI);
    prob.addMatrixOperator(opDiv, dow, i); 
  }
  
  for (int i = 0; i < dow; ++i) {
    // ===== add boundary conditions =====
    prob.addDirichletBC(1, i, i, new G(i, vel)); // inflow
    
    for (BoundaryType nr = 3; nr <= 5; ++nr)
      prob.addDirichletBC(nr, i, i, new Constant(0.0)); // wall
  }
  
  // ===== set initial solution =====
  for (int i = 0; i < dow; ++i)
    *u[i] << eval(new G(i, vel));
  
  AdaptInstationary adaptInstat("adapt", prob, adaptInfo, instat, adaptInfo);
  adaptInstat.adapt();

  AMDiS::finalize();
}
