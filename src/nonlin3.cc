#include "AMDiS.h"
#include "expressions/diff_expr.hpp"

struct u_; // create a class to identify the unknown in the expressions

using namespace AMDiS;

// solve: -laplace(u) + u^4 = f(x) in Omega,    u = g on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemNonLin prob("nonlin");
  prob.initialize(INIT_ALL);

  FiniteElemSpace const* V = prob.getFeSpace(0);
  DOFVector<double>& u_k =  *prob.getSolution(0);
  
  auto f = X()*X() + 1.0;
  auto coeff = pow<4>(valueOf<u_>(u_k)) - f;

  // J_u L(u_k)[du]
  Operator opJL(V, V);
    addSOT(opJL, 1.0);
    addZOT(opJL, diff<u_>(coeff) ); 
  prob.addMatrixOperator(opJL, 0, 0);

  // -L(u_k)
  Operator opL(V);
    addSOT(opL, -1.0);
  opL.setUhOld(&u_k);
  prob.addVectorOperator(opL, 0);

  Operator opF(V);
    addZOT(opF, -coeff);
  prob.addVectorOperator(opF, 0);

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new Constant(0.0));

  // ===== start adaption loop =====
  AdaptInfo adaptInfo("adapt");

  u_k << 0.0;
  prob.solve(&adaptInfo);
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);
  
  std::cout << coeff << "\n"; 
  std::cout << diff<u_>(coeff) << "\n"; 
  std::cout << diff<2,u_>(coeff) << " --> " << simplify(diff<2,u_>(coeff)) << "\n"; 
  std::cout << diff<3,u_>(coeff) << " --> " << simplify(diff<3,u_>(coeff)) << "\n"; 
  std::cout << diff<4,u_>(coeff) << " --> " << simplify(diff<4,u_>(coeff)) << "\n"; 

  AMDiS::finalize();
}
