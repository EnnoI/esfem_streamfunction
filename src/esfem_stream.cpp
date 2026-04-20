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
#include <amdis/esfem_stream/StreamFunctionOperator.hpp>
#include <amdis/esfem_stream/StreamFunctionRhsOperator.hpp>

#include <src/data/stokes.hpp>

auto const _phi = Dune::Indices::_0;
auto const _psi = Dune::Indices::_1;
auto const _omega = Dune::Indices::_2;
auto const _vn = Dune::Indices::_3;
auto const _pn = Dune::Indices::_4;
auto const _Hn = Dune::Indices::_5;

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

  int constexpr kg = 4; // order of geometry
  int constexpr ku = 2; // lagrange order of esfem streamfunction problem

  auto surfaceBasisFactory = power<3>(lagrange<kg>(), blockedInterleaved());
  GlobalBasis hostSurfaceBasis{adaptiveHostGridPtr->leafGridView(), surfaceBasisFactory};
  auto surfaceFctDOF = makeDOFVector(hostSurfaceBasis);
  auto surfaceFct = valueOf(surfaceFctDOF);
  surfaceFct << [](FieldVector<double,3> const& x) { return x / x.two_norm(); };

  using CurvedGrid = Dune::CurvedGrid<AdaptiveHostGrid, decltype(surfaceFct)>;
  using Grid = AdaptiveGrid<CurvedGrid>;

  Dune::CurvedGrid grid(*adaptiveHostGridPtr, surfaceFct);
  AdaptiveGrid<CurvedGrid> adaptiveGrid(grid);

  // interpolate identity on curved grid
  GlobalBasis surfaceBasis{adaptiveGrid.leafGridView(), surfaceBasisFactory};
  auto surfaceIdentityFctDOF = makeDOFVector(surfaceBasis);
  auto surfaceIdentityFct = valueOf(surfaceIdentityFctDOF);

  // Define problems
  using Prob = ProblemStat<LagrangeBasis<Grid,ku,ku,ku,ku,ku,ku>>; // phi,psi,omega,vn,pn,Hn
  Prob prob("esfem_stream", adaptiveGrid);
  prob.initialize(INIT_ALL);

  auto bgnBasisFactory = composite(surfaceBasisFactory, lagrange<kg>()); // Y, H
  ProblemStat bgnProb("bgn", adaptiveGrid, bgnBasisFactory);
  bgnProb.initialize(INIT_ALL);

  // General parameters
  double const dt = Parameters::get<double>("parameters->dt").value_or(1e-3);
  int const kh = Parameters::get<int>("parameters->kh").value_or(2);
  double surface = 1.0; // dummy argument
  double const alpha = Parameters::get<double>("parameters->f_T").value_or(1.0);
  double const c = Parameters::get<double>("parameters->c").value_or(0.95);
  double const d = Parameters::get<double>("parameters->d").value_or(0.96);

  // ---== Define streamfunction problem ==---
  std::array<double,3> data;
  auto Hgf = makeGridFunction(prob.solution(_Hn), prob.globalBasis()->gridView());
  prob.addMatrixOperator( makeOperator(tag::stream_function{kg, kh, surface, data, Hgf}, 1.0) );

  auto f0 = F{c,d};
  auto f = alpha * Dune::analyticGridFunction<Grid>(f0);
  prob.addVectorOperator( makeOperator(tag::streamfunction_rhs{kh, surface, Hgf}, f, 20) );

  // ---== Define BGN problem ==---
  // <H*n, y> + <grad Y, grad y> = -<grad X, grad y>
  // <Y * n, h> = <vn*dt, h>
  auto vn = [](FieldVector<double,3> const& x)->double { return std::cos(x[0]+x[1]); };
  bgnProb.addVectorOperator( makeOperator(tag::bgnrhs{kh, dt, surface}, prob.solution(_vn), 20) );
  bgnProb.addMatrixOperator( makeOperator(tag::bgn{kh, surface}, 1.0) );

  int order_writer = Parameters::get<int>("surface->write files->order").value_or(2);
  Dune::Vtk::LagrangeDataCollector dataCollector(bgnProb.gridView(), order_writer);
  Dune::VtkUnstructuredGridWriter writer(dataCollector);

  adaptiveGrid.globalRefine(Parameters::get<int>("postprocessing->level").value_or(1));
  // bgnProb.globalRefine(Parameters::get<int>("postprocessing->level").value_or(1));

  surfaceFct << [](FieldVector<double,3> const& x) { return x / x.two_norm(); };
  surfaceIdentityFct << [](FieldVector<double,3> const& x) { return x; };
  prob.solution(_Hn) << [](FieldVector<double,3> const& x) { return -2.0; };
  
  // Add point data to writer
  {
    writer.addPointData(prob.solution(_phi), Dune::Vtk::FieldInfo{"phi", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer.addPointData(prob.solution(_psi), Dune::Vtk::FieldInfo{"psi", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer.addPointData(prob.solution(_omega), Dune::Vtk::FieldInfo{"omega", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer.addPointData(prob.solution(_vn), Dune::Vtk::FieldInfo{"vn", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer.addPointData(prob.solution(_pn), Dune::Vtk::FieldInfo{"pn", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer.addPointData(prob.solution(_Hn), Dune::Vtk::FieldInfo{"Hn", 1, Dune::Vtk::RangeTypes::SCALAR});

    writer.addPointData(bgnProb.solution(_H), Dune::Vtk::FieldInfo{"H", 1, Dune::Vtk::RangeTypes::SCALAR});
    writer.addPointData(bgnProb.solution(_Y), Dune::Vtk::FieldInfo{"Y", 3, Dune::Vtk::RangeTypes::VECTOR});
  }

  // Assemble and solve linear systems
  AdaptInfo adaptInfo("adapt");

  // ---== Solve stream function problem ==---
  prob.assemble(adaptInfo);
  prob.solve(adaptInfo);

  // ---== Solve BGN problem ==---
  bgnProb.assemble(adaptInfo);
  bgnProb.solve(adaptInfo);

  writer.write("output/esfem1");

  // Evolve surface
  auto& X = surfaceFct.coefficients().impl().vector();
  auto dXDOF = makeDOFVector(surfaceBasis);
  auto dXFct = valueOf(dXDOF);
  dXFct << bgnProb.solution(_Y);
  auto const& dx = dXFct.coefficients().impl().vector();
  X += dx;

  writer.write("output/esfem2");

}
