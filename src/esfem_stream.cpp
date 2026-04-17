#include <amdis/AMDiS.hpp>
#include <amdis/ProblemStat.hpp>

#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/curvedgrid/gridfunctions/discretegridviewfunction.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/3d/alugrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

using namespace AMDiS;
int main(int argc, char** argv)
{
  Environment env(argc, argv);

  // using HostGrid = Dune::FoamGrid<2,3>;
  using HostGrid = Dune::ALUGrid<2,3,Dune::simplex,Dune::conforming>;
  using AdaptiveHostGrid = AdaptiveGrid<HostGrid>;
  using namespace Dune::Functions::BasisFactory;

  const std::string meshfile = "macro/sphere.msh";
  Dune::GmshReader<HostGrid> reader(meshfile); 
  std::unique_ptr<HostGrid> hostGridPtr = reader.createGrid();
  std::unique_ptr<AdaptiveHostGrid> adaptiveHostGridPtr = std::make_unique<AdaptiveHostGrid>(*hostGridPtr);

  auto surfaceBasisFactory = power<3>(lagrange<3>(), blockedInterleaved());
  GlobalBasis surfaceBasis{adaptiveHostGridPtr->leafGridView(), surfaceBasisFactory};
  auto surfaceFctDOF = makeDOFVector(surfaceBasis);
  // auto surfaceFct = valueOf(surfaceFctDOF);
  // surfaceFct << [](FieldVector<double,3> const& x) { return x / x.two_norm(); };

  auto surfaceFct = Dune::discreteGridViewFunction<3>(adaptiveHostGridPtr->leafGridView(), 3);

  using CurvedGrid = Dune::CurvedGrid<AdaptiveHostGrid, decltype(surfaceFct)>;

  Dune::CurvedGrid grid(*adaptiveHostGridPtr, surfaceFct);
  AdaptiveGrid<CurvedGrid> adaptiveGrid(grid);

  auto bgnBasisFactory = composite(power<3>(lagrange<3>(), blockedInterleaved()), lagrange<3>());
  ProblemStat bgnProb("bgn", adaptiveGrid, bgnBasisFactory);
  bgnProb.initialize(INIT_ALL);

  // bgnProb.globalRefine(1);
}
