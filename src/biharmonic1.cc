#include "AMDiS.h"

using namespace AMDiS;

// solve: -laplace(u) = f(x) in Omega,    u = g on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("biharmonic");
  prob.initialize(INIT_ALL);

  double nu = 1.0;
  Parameters::get("biharmonic->nu", nu);
  
  // ===== define operators =====
  Operator opL_phi(prob.getFeSpace(0), prob.getFeSpace(1));
  addSOT(opL_phi, nu);
  
  Operator opF(prob.getFeSpace(0));
  // F(x) = (y, -x)  =>  f(x) = -2
  addZOT(opF, -2.0);
  
  Operator opL_psi(prob.getFeSpace(1), prob.getFeSpace(0));
  addSOT(opL_psi, 1.0); 
  
  Operator opM_phi(prob.getFeSpace(1), prob.getFeSpace(1));
  addZOT(opM_phi, 1.0);
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opL_phi, 0, 1);   // -nu*laplace(phi)
  prob.addMatrixOperator(opL_psi, 1, 0);   // -laplace(psi)
  prob.addMatrixOperator(opM_phi, 1, 1);   // phi
  prob.addVectorOperator(opF, 0);          // f(x)

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 1, 0, new Constant(0.0));
  prob.addDirichletBC(nr, 0, 1, new Constant(0.0));

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt", prob.getNumComponents());
  
  // ===== assemble and solve linear system =====
  AdaptStationary adaptStat("adapt", prob, adaptInfo);
  adaptStat.adapt();
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);

  AMDiS::finalize();
}
