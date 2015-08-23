#include "AMDiS.h"

using namespace AMDiS;

struct G : AbstractFunction<double, WorldVector<double> >
{
  double operator()(WorldVector<double> const& x) const 
  {
    return std::exp(x*x) + std::cos(10.0 * (x[0] + x[1]));
  }
};

// solve: dt(u) - laplace(u) = f(x) in Omega,    u = 0 on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("heat");
  prob.initialize(INIT_ALL);
  
  ProblemInstat instat("heat", prob);
  instat.initialize(INIT_NOTHING);

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt");

  // ===== define operators =====
  Operator opMatrix(prob.getFeSpace(), prob.getFeSpace());
  addZOT(opMatrix, 1.0/ref_(adaptInfo.getTimestepPtr())); // 1/tau * u
  addSOT(opMatrix, 1.0);  // -laplace(u)
  
  DOFVector<double>& u = *prob.getSolution(0);
  
  Operator opVector(prob.getFeSpace());
  addZOT(opVector, valueOf(u) / ref_(adaptInfo.getTimestepPtr())); // 1/tau * u_old
  addZOT(opVector, dot(X(), X()) + 10.0); // f(x) = x*x + 10
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opMatrix, 0, 0);   
  prob.addVectorOperator(opVector, 0);

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new G); // g(x) = exp(x*x) + cos(10*(x_0 + x_1)))
  
  // ===== set initial solution =====
  u << constant(0.0); // u(x,t=0) := 0
  
  AdaptInstationary adaptInstat("adapt", prob, adaptInfo, instat, adaptInfo);
  adaptInstat.adapt();

  AMDiS::finalize();
}
