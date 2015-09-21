#include "AMDiS.h"

using namespace AMDiS;

// solve: -laplace(u) + u^4 = f(x) in Omega,    u = 0 on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("nonlin");
  prob.initialize(INIT_ALL);

  FiniteElemSpace const* V = prob.getFeSpace(0);
  DOFVector<double>& u_k =  *prob.getSolution(0);

  // J_u L(u_k)[du]
  Operator opJL(V, V);
    addSOT(opJL, 1.0);
    addZOT(opJL, 4.0 * pow<3>(valueOf(u_k)) ); 
  prob.addMatrixOperator(opJL, 0, 0);

  // -L(u_k)
  Operator opL(V);
    addSOT(opL, -1.0);
  opL.setUhOld(&u_k);
  prob.addVectorOperator(opL, 0);

  Operator opF(V);
    addZOT(opF, -pow<4>(valueOf(u_k)) + X()*X() + 1.0 );
  prob.addVectorOperator(opF, 0);

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new Constant(0.0));

  // ===== start adaption loop =====
  AdaptInfo adaptInfo("adapt");
  AdaptStationary adaptStat("adapt", prob, adaptInfo);

  u_k << 0.0;
  for (int k = 0; k < adaptInfo.getMaxSpaceIteration(); ++k) {
    adaptInfo.reset();
    adaptStat.adapt();
  }
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);

  AMDiS::finalize();
}
