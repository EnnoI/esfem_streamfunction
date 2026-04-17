#include <amdis/AMDiS.hpp>

#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

using namespace AMDiS;
int main(int argc, char** argv)
{
  Environment env(argc, argv);

  using HostGrid = Dune::FoamGrid<2,3>;
  using AdaptiveHostGrid = AdaptiveGrid<HostGrid>;
  using namespace Dune::Functions::BasisFactory;

  const std::string meshfile = "macro/sphere.msh";
  Dune::GmshReader<HostGrid> reader(meshfile); 
  std::unique_ptr<HostGrid> hostGridPtr = reader.createGrid();
  std::unique_ptr<AdaptiveHostGrid> adaptiveHostGridPtr = std::make_unique<AdaptiveHostGrid>(*hostGridPtr);

  auto surfaceBasisFactory = power<3>(lagrange<3>(), blockedInterleaved());
  GlobalBasis surfaceBasis{adaptiveHostGridPtr->leafGridView(), surfaceBasisFactory};
  auto surfaceFctDOF = makeDOFVector(surfaceBasis);
  auto surfaceFct = valueOf(surfaceFctDOF);
  surfaceFct << [](FieldVector<double,3> const& x) { return x; };


  Dune::CurvedGrid(*adaptiveHostGridPtr, surfaceFct);

  // your code comes here
}
