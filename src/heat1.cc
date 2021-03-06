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

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt");

  // ===== define operators =====
  Operator opMatrix(prob.getFeSpace(), prob.getFeSpace());
  addZOT(opMatrix, 1.0/adaptInfo.getTimestep()); // fixed timestep
  addSOT(opMatrix, 1.0);
  
  DOFVector<double>& u = *prob.getSolution(0);
  
  Operator opVector(prob.getFeSpace());
  addZOT(opVector, valueOf(u) / adaptInfo.getTimestep());
  addZOT(opVector, dot(X(), X()) + 10.0); // f(x) = x*x + 10
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opMatrix, 0, 0);   // -laplace(u)
  prob.addVectorOperator(opVector, 0);            // f(x)

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new G); // g(x) = exp(x*x) + cos(10*(x_0 + x_1)))
  
  // ===== set initial solution =====
  u << constant(0.0); // u(x,t=0) := 0
  
  while (adaptInfo.getTime() < adaptInfo.getEndTime()) {
    adaptInfo.setTime(adaptInfo.getTime() + adaptInfo.getTimestep()); // t += tau
    
    // ===== assemble and solve linear system =====
    prob.assemble(&adaptInfo);
    prob.solve(&adaptInfo);
    
    // ===== write solution to file =======
    prob.writeFiles(&adaptInfo, false);
  }

  AMDiS::finalize();
}
