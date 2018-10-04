// Exercise3: system of equations equation
//
// Reference: see https://goo.gl/JK8EUI for a detailed description of possible 
//            expression terms.
//
// Compile-and-run: 
// > cd build
// > make exercise3
// > cd ..
// > build/exercise3 init/exercise3.dat.2d
//

#include "AMDiS.h"

using namespace AMDiS;

// solve: -laplace(u) = f(x) in Omega,    u = g on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("biharmonic");
  prob.initialize(INIT_ALL);
  
  // ===== define operators =====
  Operator opL_v(prob.getFeSpace(0), prob.getFeSpace(1));
  addSOT(opL_v, -1.0);
  
  Operator opM_u(prob.getFeSpace(0), prob.getFeSpace(0));
  addZOT(opM_u, 1.0);
  
  Operator opF(prob.getFeSpace(0));
  // F(x) = (y, -x)  =>  f(x) = -2
  addZOT(opF, eval(new Random(0.0, 1.0)) );
  
  double epsilon = 0.1;
  Parameters::get("epsilon", epsilon);
  
  Operator opL_u(prob.getFeSpace(1), prob.getFeSpace(0));
  addSOT(opL_u, -epsilon); 
  addZOT(opL_u,  1.0);
  
  Operator opM_v(prob.getFeSpace(1), prob.getFeSpace(1));
  addZOT(opM_v, 1.0);
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opM_u, 0, 0);   // phi
  prob.addMatrixOperator(opL_v, 0, 1);   // -nu*laplace(phi)
  prob.addMatrixOperator(opL_u, 1, 0);   // -laplace(psi)
  prob.addMatrixOperator(opM_v, 1, 1);   // phi
  prob.addVectorOperator(opF, 0);        // f(x)

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt", prob.getNumComponents());
  
  // ===== assemble and solve linear system =====
  AdaptStationary adaptStat("adapt", prob, adaptInfo);
  adaptStat.adapt();
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);

  AMDiS::finalize();
}
