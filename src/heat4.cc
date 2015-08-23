#include "AMDiS.h"
#include "BaseProblem.h"

using namespace AMDiS;

struct G : AbstractFunction<double, WorldVector<double> >
{
  double operator()(WorldVector<double> const& x) const 
  {
    return std::exp(x*x) + std::cos(10.0 * (x[0] + x[1]));
  }
};

struct Heat : StandardBaseProblem
{
  Heat(std::string name) 
    : StandardBaseProblem(name),
      u(NULL)
  {}
  
  void initData()
  {
    u = getProblem()->getSolution(0);
  }
  
  // ===== set initial solution =====
  void solveInitialProblem(AdaptInfo* adaptInfo)
  {
    *u << constant(0.0); // u(x,t=0) := 0
  }
  
  void closeTimestep(AdaptInfo* adaptInfo)
  { FUNCNAME("Heat::closeTimestep");
    
    MSG("int{ u } = %e\n", integrate( valueOf(u) ));
  }
  
  // ===== define operators =====
  void fillOperators()
  {
    Operator* opMatrix = new Operator(getFeSpace(), getFeSpace());
    addSOT(opMatrix, 1.0);  // -laplace(u)
    
    Operator* opVector = new Operator(getFeSpace());
    addZOT(opVector, dot(X(), X()) + 10.0); // f(x) = x*x + 10
    
    // ===== add operators to problem =====
    addTimeOperator(getProblem(), 0, 0); // 1/tau*u, 1/tau*u_old
    getProblem()->addMatrixOperator(opMatrix, 0, 0);   
    getProblem()->addVectorOperator(opVector, 0);
  }
  
  // ===== add boundary conditions =====
  void fillBoundaryConditions()
  {
    BoundaryType nr = 1;
    getProblem()->addDirichletBC(nr, 0, 0, new G); // g(x) = exp(x*x) + cos(10*(x_0 + x_1)))
  }
  
private:
  DOFVector<double>* u;
};

// solve: dt(u) - laplace(u) = f(x) in Omega,    u = 0 on Gamma
int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  Heat prob("heat");
  prob.initialize(INIT_ALL);
  prob.initTimeInterface(); // calls initData(), fillOperators(), ...

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt");
  
  AdaptInstationary adaptInstat("adapt", prob, adaptInfo, prob, adaptInfo);
  adaptInstat.adapt();

  AMDiS::finalize();
}
