#include "AMDiS.h"

using namespace AMDiS;

// solve: -laplace(u) = f(x) in Omega,    u = 0 on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("poisson");
  prob.initialize(INIT_ALL);

  // ===== define operators =====
  Operator opLaplace(prob.getFeSpace(), prob.getFeSpace());
  addSOT(opLaplace, 1.0); // <grad(u), grad(theta)>
  
  Operator opF(prob.getFeSpace());
  addZOT(opF, 1.0); // <f(x), theta>
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opLaplace, 0, 0);
  prob.addVectorOperator(opF, 0);

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new Constant(0.0));

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt");
  
  // ===== assemble and solve linear system =====
  prob.assemble(&adaptInfo);
  prob.solve(&adaptInfo);
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);
  
  AMDiS::finalize();
}
