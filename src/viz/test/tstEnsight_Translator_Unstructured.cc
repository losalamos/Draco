//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/tstEnsight_Translator_Unstructured.cc
 * \author Thomas M. Evans, Ryan T. Wollaeger
 * \date   Wednesday, Oct 03, 2018, 15:27 pm
 * \brief  Ensight_Translator unstructured mesh test.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/path.hh"
#include "viz/Ensight_Translator.hh"

using rtt_viz::Ensight_Translator;

template <typename IT>
void ensight_dump_test_unstr2d(rtt_dsxx::UnitTest &ut, bool const binary) {

  // short-cuts
  typedef std::vector<std::string> vec_s;
  typedef std::vector<IT> vec_i;
  typedef std::vector<vec_i> vec2_i;
  typedef std::vector<double> vec_d;
  typedef std::vector<vec_d> vec2_d;

  if (binary)
    std::cout << "\nGenerating binary files...\n" << std::endl;
  else
    std::cout << "\nGenerating ascii files...\n" << std::endl;

  // >>> SET SCALAR ENSIGHT CTOR ARGS

  std::string prefix = "unstr2d_testproblem";
  if (binary)
    prefix += "_binary";

  int icycle = 1;
  double time = .01;
  double dt = .01;
  const bool static_geom = false;

  std::string const gd_wpath = rtt_dsxx::getFilenameComponent(
      ut.getTestInputPath(), rtt_dsxx::FC_NATIVE);

  // >>> INITIALIZE AND SET VECTOR DATA

  // dimensions
  size_t ncells = 2;
  size_t nvert = 8;
  size_t ndim = 2;
  size_t ndata = 2;
  std::vector<size_t> nvert_per_cell = {4, 6};

  // size the cell-vertex vector
  vec2_i ipar = {vec_i(nvert_per_cell[0]), vec_i(nvert_per_cell[1])};
  Check(ipar.size() == ncells);
  vec2_d vrtx_data(nvert, vec_d(ndata));
  vec2_d cell_data(ncells, vec_d(ndata));
  vec2_d pt_coor(nvert, vec_d(ndim));
  // set the element type to be unstructured
  vec_i iel_type(ncells, rtt_viz::unstructured);
  // set the vertex-centered and cell-centered data names
  vec_s vdata_names = {"Densities", "Temperatures"};
  Check(vdata_names.size() == ndata);
  vec_s cdata_names = {"Velocity", "Pressure"};
  Check(cdata_names.size() == ndata);
  // region information
  vec_i rgn_index = {1, 1};
  Check(rgn_index.size() == ncells);
  vec_s rgn_name = {"RGN_A"};
  vec_i rgn_data = {1};

  // create some arbitrary cell and vertex based data (as in structured test)
  for (size_t i = 0; i < ndata; i++) {
    // cell data
    for (size_t cell = 0; cell < ncells; cell++) {
      Check(1 + cell < INT_MAX);
      cell_data[cell][i] = static_cast<int>(1 + cell);
    }

    // vrtx data
    for (size_t v = 0; v < nvert; v++) {
      Check(1 + v < INT_MAX);
      vrtx_data[v][i] = static_cast<int>(1 + v);
    }
  }

  // Build path for the input file "cell_data_unstr2d"
  std::string const cdInputFile =
      ut.getTestSourcePath() + std::string("cell_data_unstr2d");
  std::ifstream input(cdInputFile.c_str());
  if (!input)
    ITFAILS;

  for (size_t i = 0; i < pt_coor.size(); i++) {
    for (size_t j = 0; j < pt_coor[i].size(); j++)
      input >> pt_coor[i][j];
  }
  for (size_t i = 0; i < ipar.size(); i++) {
    for (size_t j = 0; j < ipar[i].size(); j++)
      input >> ipar[i][j];
  }

  // build an Ensight_Translator (make sure it overwrites any existing stuff)
  Ensight_Translator translator(prefix, gd_wpath, vdata_names, cdata_names,
                                true, static_geom, binary);

  translator.ensight_dump(icycle, time, dt, ipar, iel_type, rgn_index, pt_coor,
                          vrtx_data, cell_data, rgn_data, rgn_name);

  std::vector<double> dump_times = translator.get_dump_times();
  if (dump_times.size() != 1)
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(dump_times[0], 0.01))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("ensight_dump_test_unstr2d finished successfully.");
  else
    FAILMSG("ensight_dump_test_unstr2d did not finish successfully.");
  return;
}

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // ASCII dumps
    bool binary(false);
    ensight_dump_test_unstr2d<int>(ut, binary);

    // Binary dumps
    binary = true;
    ensight_dump_test_unstr2d<int>(ut, binary);

    // ASCII dumps with unsigned integer data
    binary = false;
    ensight_dump_test_unstr2d<uint32_t>(ut, binary);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of viz/test/tstEnsight_Translator_Unstructured.cc
//----------------------------------------------------------------------------//
