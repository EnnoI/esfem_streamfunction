#include <string>
#include <cmath>

#include <amdis/AMDiS.hpp>
#include <amdis/ProblemStat.hpp>
#include <amdis/AdaptiveGrid.hpp>
#include <amdis/DOFVector.hpp>

#include <dune/curvedgrid/curvedgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/alugrid/3d/alugrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/vtk/vtkreader.hh>
#include <dune/vtk/gridcreators/lagrangegridcreator.hh>
#include <dune/vtk/vtkwriter.hh>
#include <dune/vtk/datacollectors/lagrangedatacollector.hh>

#include <amdis/esfem_stream/BgnOperator.hpp>
#include <amdis/esfem_stream/BgnRhsOperator.hpp>

auto const _Y = Dune::Indices::_0;
auto const _H = Dune::Indices::_1;

using namespace AMDiS;
int main(int argc, char** argv)
{
  Environment env(argc, argv);

  // using HostGrid = Dune::FoamGrid<2,3>;
  using HostGrid = Dune::ALUGrid<2,3,Dune::simplex,Dune::conforming>;
  using AdaptiveHostGrid = AdaptiveGrid<HostGrid>;
  using namespace Dune::Functions::BasisFactory;

  const std::string meshfile = "macro/sphere.msh";
  using Reader = Dune::GmshReader<HostGrid>;
  auto hostGridPtr = Reader::read(meshfile,true,false);

  // const std::string meshfile = "macro/h025a.vtu";
  // std::unique_ptr<HostGrid> hostGridPtr = Dune::Vtk::VtkReader<HostGrid, Dune::Vtk::LagrangeGridCreator<HostGrid>>::createGridFromFile(meshfile);

  std::unique_ptr<AdaptiveHostGrid> adaptiveHostGridPtr = std::make_unique<AdaptiveHostGrid>(*hostGridPtr);

  int constexpr kg = 4;

  auto surfaceBasisFactory = power<3>(lagrange<kg>(), blockedInterleaved());
  GlobalBasis surfaceBasis{adaptiveHostGridPtr->leafGridView(), surfaceBasisFactory};
  auto surfaceFctDOF = makeDOFVector(surfaceBasis);
  auto surfaceFct = valueOf(surfaceFctDOF);
  surfaceFct << [](FieldVector<double,3> const& x) { return x / x.two_norm(); };

  using CurvedGrid = Dune::CurvedGrid<AdaptiveHostGrid, decltype(surfaceFct)>;

  Dune::CurvedGrid grid(*adaptiveHostGridPtr, surfaceFct);
  AdaptiveGrid<CurvedGrid> adaptiveGrid(grid);

  // interpolate identity on curved grid
  GlobalBasis surfaceIdentityBasis{adaptiveGrid.leafGridView(), surfaceBasisFactory};
  auto surfaceIdentityFctDOF = makeDOFVector(surfaceIdentityBasis);
  auto surfaceIdentityFct = valueOf(surfaceIdentityFctDOF);

  auto bgnBasisFactory = composite(surfaceBasisFactory, lagrange<kg>());
  ProblemStat bgnProb("bgn", adaptiveGrid, bgnBasisFactory);
  bgnProb.initialize(INIT_ALL);

  double dt = Parameters::get<double>("parameters->dt").value_or(1e-3);
  int kh = Parameters::get<int>("parameters->kh").value_or(2);
  double surface = 1.0;

  // <H*n, y> + <grad Y, grad y> = -<grad X, grad y>
  // <Y * n, h> = <vn*dt, h>
  auto vn = [](FieldVector<double,3> const& x)->double { return std::cos(x[0]+x[1]); };
  // bgnProb.addVectorOperator( bgnrhs(kh, dt, surface, vn, gradientOf(surfaceFct)) ); // _cache version
  bgnProb.addVectorOperator( makeOperator(tag::bgnrhs{kh, dt, surface}, vn, 20) );
  bgnProb.addMatrixOperator( makeOperator(tag::bgn{kh, surface}, 1.0) );

  int order_writer = Parameters::get<int>("surface->write files->order").value_or(2);
  Dune::Vtk::LagrangeDataCollector dataCollector(bgnProb.gridView(), order_writer);
  Dune::VtkUnstructuredGridWriter writer(dataCollector);

  bgnProb.globalRefine(Parameters::get<int>("postprocessing->level").value_or(1));

  surfaceFct << [](FieldVector<double,3> const& x) { return x / x.two_norm(); };
  surfaceIdentityFct << [](FieldVector<double,3> const& x) { return x; };

  writer.addPointData(bgnProb.solution(_H), Dune::Vtk::FieldInfo{"H", 1, Dune::Vtk::RangeTypes::SCALAR});
  writer.addPointData(bgnProb.solution(_Y), Dune::Vtk::FieldInfo{"Y", 3, Dune::Vtk::RangeTypes::VECTOR});

  // Assemble and solve linear systems
  AdaptInfo adaptInfo("adapt");

  bgnProb.assemble(adaptInfo);
  bgnProb.solve(adaptInfo);

  writer.write("output/esfem1");

  auto& X = surfaceFct.coefficients().impl().vector();
  auto const& dx = bgnProb.solution(_Y).coefficients().impl().vector();
  X += dx;

  writer.write("output/esfem2");

}
