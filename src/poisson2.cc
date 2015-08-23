#include "AMDiS.h"

using namespace AMDiS;

struct G : AbstractFunction<double, WorldVector<double> >
{
  double operator()(WorldVector<double> const& x) const 
  {
    return std::exp(x*x) + std::cos(10.0 * (x[0] + x[1]));
  }
};

// solve: -laplace(u) = f(x) in Omega,    u = g on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("poisson");
  prob.initialize(INIT_ALL);

  // ===== define operators =====
  Operator opLaplace(prob.getFeSpace(), prob.getFeSpace());
  addSOT(opLaplace, 1.0);
  
  Operator opF(prob.getFeSpace());
  addZOT(opF, dot(X(), X()) + 1.0); // f(x) = x*x + 1
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opLaplace, 0, 0);   // -laplace(u)
  prob.addVectorOperator(opF, 0);            // f(x)

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new G); // g(x) = exp(x*x) + cos(10(x_0 + x_1)))

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt");
  
  // ===== assemble and solve linear system =====
  prob.assemble(&adaptInfo);
  prob.solve(&adaptInfo);
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);

  AMDiS::finalize();
}
