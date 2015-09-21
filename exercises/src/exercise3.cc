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

int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("ch");
  prob.initialize(INIT_ALL);
  
  // ===== define operators =====
  
  // operators for equation 0
  Operator opL_v(prob.getFeSpace(0), prob.getFeSpace(1));
  Operator opM_u(prob.getFeSpace(0), prob.getFeSpace(0));
  Operator opF(prob.getFeSpace(0));
  
  // add operator terms to the operators
  
  // operators for equation 1
  Operator opL_u(prob.getFeSpace(1), prob.getFeSpace(0));
  Operator opM_v(prob.getFeSpace(1), prob.getFeSpace(1));
  
  // add operator terms to the operators
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opM_u, 0, 0);   // u
  prob.addMatrixOperator(opL_v, 0, 1);   // laplace(v)
  prob.addMatrixOperator(opL_u, 1, 0);   // laplace(u)
  prob.addMatrixOperator(opM_v, 1, 1);   // v
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
