// Example: exercise2-advanced
//
// Goal: implement a simple adaptive strategy based on exact solution
//
// Reference: see https://fusionforge.zih.tu-dresden.de/plugins/mediawiki/wiki/amdis/index.php/AMDiS::Expressions
//            for a detailed description of possible expression terms.
//
// Compile-and-run: 
// > cd build
// > make exercise2b
// > cd ..
// > build/exercise2b init/exercise2.dat.2d
//

#include "AMDiS.h"

#include <array>
#include <boost/numeric/mtl/mtl.hpp>

using namespace AMDiS;


void convergence_factor(std::vector<std::array<double,2>> const& error)
{
  mtl::dense_vector<double> rhs_res(error.size());
  mtl::dense2D<double> A_res(error.size(), 2);
  for (size_t i = 0; i < error.size(); ++i) {
    rhs_res(i) = std::log(error[i][1]);
    A_res(i, 0) = 1;
    A_res(i, 1) = std::log(error[i][0]);
  }
  
  mtl::dense_vector<double> convergence(2);
  
  mtl::dense2D<double> Q(num_rows(A_res), num_rows(A_res)), R(num_rows(A_res), num_cols(A_res));
  boost::tie(Q, R)= mtl::matrix::qr(A_res);
  
  mtl::dense_vector<double> b(trans(Q)*rhs_res);
  mtl::irange rows(0,num_cols(A_res));
  convergence = mtl::matrix::upper_trisolve(R[rows][rows], b);
  
  std::cout << "\n|u - u_h| = C * h^k\n   C = " << std::exp(convergence(0)) << "\n   k = " << convergence(1) << "\n\n";
}

inline double calcMeshSizes(Mesh* mesh) 
{
  TraverseStack stack;
  ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
  double maxH = 0.0;
  while (elInfo) {
    auto coords = elInfo->getCoords();
    for (int i = 0; i < coords.getSize(); i++)
      for (int j = i+1; j < coords.getSize(); j++)
	  maxH = std::max(maxH, norm(coords[i] - coords[j]));
    elInfo = stack.traverseNext(elInfo);
  }
  return maxH;
}

struct G : AbstractFunction<double, WorldVector<double> >
{
  double operator()(WorldVector<double> const& x) const 
  {
    return std::exp(x*x) + std::cos(10.0 * (x[0] + x[1]));
  }
};

// solve: -laplace(u) = f(x) in Omega,    u = g on Gamma
int main(int argc, char* argv[])
{ FUNCNAME("Main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat prob("poisson");
  prob.initialize(INIT_ALL);

  // ===== define operators =====
  Operator opLaplace(prob.getFeSpace(), prob.getFeSpace());
  addSOT(opLaplace, 1.0);
  
  Operator opF(prob.getFeSpace());
  auto f = dot(X(), X()) + 1.0;
  addZOT(opF, f); // f(x)
  
  // ===== add operators to problem =====
  prob.addMatrixOperator(opLaplace, 0, 0);   // -laplace(u)
  prob.addVectorOperator(opF, 0);            // f(x)

  // ===== add boundary conditions =====
  BoundaryType nr = 1;
  prob.addDirichletBC(nr, 0, 0, new G); // g(x)

  // ===== create info-object, that holds parameters ===
  AdaptInfo adaptInfo("adapt");
  
  DOFVector<double>& U = *prob.getSolution(0);
  DOFVector<double> ErrVec(U), UExact(U);
  
  auto u_exact = exp(X()*X()) + cos(10.0*(X(0) + X(1)));
  
  std::vector<std::array<double, 2>> error_vec_L2;
  std::vector<std::array<double, 2>> error_vec_H1;
  std::vector<std::array<double, 2>> error_vec_P;
  
  double error = 1.e10;  
  for (int i = 0; i < adaptInfo.getMaxSpaceIteration(); ++i)
  {
    // ===== assemble and solve linear system =====
    prob.assemble(&adaptInfo);
    prob.solve(&adaptInfo);
  
    ErrVec << absolute(valueOf(U) - u_exact);
    io::writeFile(ErrVec, "error_" + std::to_string(i) + ".vtu");
    
    UExact << u_exact;
    double errorL2 = std::sqrt( integrate( pow<2>(valueOf(U) - valueOf(UExact)) ) );
    double errorH1 = std::sqrt( integrate( pow<2>(valueOf(U) - valueOf(UExact)) + unary_dot(gradientOf(U) - gradientOf(UExact)) ) );
    
    double h_max = calcMeshSizes(prob.getMesh());
    
    WorldVector<double> p; p[0] = 0.5; p[1] = 0.5;    
    double errorP = std::abs(U(p) - exp(-10.0*(p*p)));
    MSG("h(%d) = %f\n", i, h_max);
    MSG("errorL2(%d) = %e\n", i, errorL2);
    MSG("errorH1(%d) = %e\n", i, errorH1);
    MSG("errorP(%d) = %e\n", i, errorP);
        
    error_vec_L2.push_back({h_max, errorL2});
    error_vec_H1.push_back({h_max, errorH1});
    error_vec_P.push_back({h_max, errorP});
    
    error = errorL2;
    if (error < adaptInfo.getSpaceTolerance(0))
      break;
    
    // refine mesh
    RefinementManager* refManager = prob.getRefinementManager();
    Flag f = refManager->globalRefine(prob.getMesh(), 1);
  }
  
  std::cout << "L2-norm:\n";
  convergence_factor(error_vec_L2);
  
  std::cout << "H1-norm:\n";
  convergence_factor(error_vec_H1);
  
  std::cout << "pointwise-norm:\n";
  convergence_factor(error_vec_P);
  
  AMDiS::finalize();
}
