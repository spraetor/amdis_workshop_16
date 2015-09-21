// Exercise2: poisson equation and postprocessing
//
// Goal: calculate norms and errors
//
// Reference: see https://goo.gl/JK8EUI for a detailed description of possible 
//            expression terms.
//
// Compile-and-run: 
// > cd build
// > make exercise2
// > cd ..
// > build/exercise2 init/exercise2.dat.2d
//

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
    
  // create some DOFVectors to store the solutions
  DOFVector<double>& U = *prob.getSolution(0);
  DOFVector<double> UExact(U);
  DOFVector<double> ERR(U);
  
  // the exact solution
  auto g = exp(X()*X()) + cos(10.0*(X(0) + X(1)));
  
  // ===== assemble and solve linear system =====
  for (int i = 0; i < 2; ++i) {
    prob.assemble(&adaptInfo);
    prob.solve(&adaptInfo);
    
    // a) Interpolate the error |uh - g| to a DOFVector
    UExact << g;
    ERR << absolute(valueOf(U) - valueOf(UExact));
    
    // b) Write the error to a file
    io::writeFile(ERR, "ERR.vtu");
    
    // c) Calculate the error-norms |uh-g|_* with * in {L1, H1}
    double errorL2 = std::sqrt( integrate( pow<2>(valueOf(U) - valueOf(UExact)) ) );
    double errorH1 = std::sqrt( integrate( pow<2>(valueOf(U) - valueOf(UExact)) + unary_dot(gradientOf(U) - gradientOf(UExact)) ) );
    std::cout << i << ") errorL2 = " << errorL2 << "\n";
    std::cout << i << ") errorH1 = " << errorH1 << "\n";
    
    // d) Refine the mesh, assemble and solve again and calculate the errors.
    RefinementManager* refManager = prob.getRefinementManager();
    Flag f = refManager->globalRefine(prob.getMesh(), 1);
  }
  
  // ===== write solution to file =======
  prob.writeFiles(&adaptInfo, true);

  AMDiS::finalize();
}
