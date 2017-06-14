//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/tstEnsightTranslator.cc
 * \author Thomas M. Evans
 * \date   Mon Jan 24 11:12:59 2000
 * \brief  Ensight_Translator test.
 * \note   Copyright (C) 2000-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/path.hh"
#include "viz/Ensight_Translator.hh"

using namespace std;
using rtt_viz::Ensight_Translator;

//---------------------------------------------------------------------------//
void ensight_dump_test(rtt_dsxx::UnitTest &ut, bool const binary) {
  if (binary)
    cout << "\nGenerating binary files...\n" << endl;
  else
    cout << "\nGenerating ascii files...\n" << endl;

  // dimensions
  int ncells = 27;
  int nvert = 64;
  int ndim = 3;
  int ndata = 2;
  int nhexvert = 8;
  int nrgn = 2;

  typedef vector<string> vec_s;
  typedef vector<int> vec_i;
  typedef vector<vec_i> vec2_i;
  typedef vector<vec2_i> vec3_i;
  typedef vector<double> vec_d;
  typedef vector<vec_d> vec2_d;
  typedef vector<vec2_d> vec3_d;

  // do an Ensight Dump
  vec2_i ipar(ncells, vec_i(nhexvert));
  vec2_d vrtx_data(nvert, vec_d(ndata, 5.0));
  vec2_d cell_data(ncells, vec_d(ndata, 10.));
  vec2_d pt_coor(nvert, vec_d(ndim));

  vec_i iel_type(ncells, rtt_viz::eight_node_hexahedron);
  vec_i rgn_index(ncells, 1);
  vec_s vdata_names(ndata, "Temperatures");
  vec_s cdata_names(ndata, "Velocity");
  vec_s rgn_name(nrgn, "RGN_A");
  vec_i rgn_data(nrgn, 1);

  // set region stuff
  rgn_name[1] = "RGN_B";
  rgn_data[1] = 2;
  for (int i = 1; i < 5; i++)
    rgn_index[i] = 2;
  rgn_index[14] = 2;
  rgn_index[15] = 2;
  rgn_index[21] = 2;
  vdata_names[1] = "Densities";
  cdata_names[1] = "Pressure";

  string prefix = "testproblem";
  if (binary)
    prefix += "_binary";

  int icycle = 1;
  double time = .01;
  double dt = .01;

  string const gd_wpath = rtt_dsxx::getFilenameComponent(ut.getTestInputPath(),
                                                         rtt_dsxx::FC_NATIVE);

  // make data
  for (int i = 0; i < ndata; i++) {
    // cell data
    for (int cell = 0; cell < ncells; cell++)
      cell_data[cell][i] = 1 + cell;

    // vrtx data
    for (int v = 0; v < nvert; v++)
      vrtx_data[v][i] = 1 + v;
  }

  // read cell data

  // Build path for the input file "cell_data"
  string const cdInputFile = ut.getTestSourcePath() + std::string("cell_data");
  ifstream input(cdInputFile.c_str());
  if (!input)
    ITFAILS;

  for (size_t i = 0; i < pt_coor.size(); i++)
    for (size_t j = 0; j < pt_coor[i].size(); j++)
      input >> pt_coor[i][j];
  for (size_t i = 0; i < ipar.size(); i++)
    for (size_t j = 0; j < ipar[i].size(); j++)
      input >> ipar[i][j];

  const bool static_geom = false;

  // Find global indices for write_part() version.

  vec2_i g_cell_indices(nrgn);
  vector<set<int>> tmp_vrtx(nrgn);
  for (int i = 0; i < ncells; i++) {
    int ipart = rgn_index[i] - 1;
    g_cell_indices[ipart].push_back(i);
    for (size_t j = 0; j < ipar[i].size(); j++)
      tmp_vrtx[ipart].insert(ipar[i][j] - 1);
  }

  typedef set<int>::const_iterator set_iter;
  vec2_i g_vrtx_indices(nrgn);
  for (int i = 0; i < nrgn; i++) {
    for (set_iter s = tmp_vrtx[i].begin(); s != tmp_vrtx[i].end(); ++s)
      g_vrtx_indices[i].push_back(*s);
  }

  // Create the equivalent data arrays for the write_part() versions.

  vec3_i p_ipar(nrgn);
  vec3_d p_vrtx_data(nrgn);
  vec3_d p_cell_data(nrgn);
  vec3_d p_pt_coor(nrgn);
  vec2_i p_iel_type(nrgn);

  for (int i = 0; i < nrgn; i++) {
    int p_ncells = g_cell_indices[i].size();
    int p_nvert = g_vrtx_indices[i].size();
    p_ipar[i].resize(p_ncells, vec_i(nhexvert));
    p_vrtx_data[i].resize(p_nvert, vec_d(ndata, 5.0));
    p_cell_data[i].resize(p_ncells, vec_d(ndata, 10.));
    p_pt_coor[i].resize(p_nvert, vec_d(ndim));
    p_iel_type[i].resize(p_ncells, rtt_viz::eight_node_hexahedron);

    for (int j = 0; j < p_nvert; j++) {
      int g = g_vrtx_indices[i][j];
      // cout << g << endl;
      p_vrtx_data[i][j] = vrtx_data[g];
      p_pt_coor[i][j] = pt_coor[g];
    }

    for (int j = 0; j < p_ncells; j++) {
      int g = g_cell_indices[i][j];
      p_cell_data[i][j] = cell_data[g];
      p_iel_type[i][j] = iel_type[g];

      for (size_t k = 0; k < ipar[g].size(); k++) {
        int tmp = ipar[g][k] - 1;

        vector<int>::iterator f =
            find(g_vrtx_indices[i].begin(), g_vrtx_indices[i].end(), tmp);

        Require(f != g_vrtx_indices[i].end());
        p_ipar[i][j][k] = f - g_vrtx_indices[i].begin() + 1;
      }
    }
  }

  // build an Ensight_Translator (make sure it overwrites any existing stuff)
  Ensight_Translator translator(prefix, gd_wpath, vdata_names, cdata_names,
                                true, static_geom, binary);

  translator.ensight_dump(icycle, time, dt, ipar, iel_type, rgn_index, pt_coor,
                          vrtx_data, cell_data, rgn_data, rgn_name);

  vec_d dump_times = translator.get_dump_times();
  if (dump_times.size() != 1)
    ITFAILS;
  if (!rtt_dsxx::soft_equiv(dump_times[0], 0.01))
    ITFAILS;

  // build another ensight translator; this should overwrite the existing
  // directories
  Ensight_Translator translator2(prefix, gd_wpath, vdata_names, cdata_names,
                                 true, static_geom, binary);

  translator2.ensight_dump(icycle, time, dt, ipar, iel_type, rgn_index, pt_coor,
                           vrtx_data, cell_data, rgn_data, rgn_name);

  // build another ensight translator from the existing dump times list;
  // thus we will not overwrite the existing directories

  Ensight_Translator translator3(prefix, gd_wpath, vdata_names, cdata_names,
                                 false, static_geom, binary);

  // now add another dump to the existing data
  translator3.ensight_dump(2, .05, dt, ipar, iel_type, rgn_index, pt_coor,
                           vrtx_data, cell_data, rgn_data, rgn_name);

  // make yet a fourth translator that will append
  Ensight_Translator translator4(prefix, gd_wpath, vdata_names, cdata_names,
                                 false, static_geom, binary);

  // add yet another dump to the existing data
  translator4.ensight_dump(3, .10, dt, ipar, iel_type, rgn_index, pt_coor,
                           vrtx_data, cell_data, rgn_data, rgn_name);

  // build an Ensight_Translator and do the per-part dump.
  string p_prefix = "part_" + prefix;
  Ensight_Translator translator5(p_prefix, gd_wpath, vdata_names, cdata_names,
                                 true, static_geom, binary);

  translator5.open(icycle, time, dt);

  for (int i = 0; i < nrgn; i++)
    translator5.write_part(i + 1, rgn_name[i], p_ipar[i], p_iel_type[i],
                           p_pt_coor[i], p_vrtx_data[i], p_cell_data[i],
                           g_vrtx_indices[i], g_cell_indices[i]);

  translator5.close();
  if (ut.numFails == 0)
    PASSMSG("ensight_dump_test finished successfully.");
  else
    FAILMSG("ensight_dump_test did not finish successfully.");
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // ASCII dumps
    bool binary(false);
    ensight_dump_test(ut, binary);

    // Binary dumps
    binary = true;
    ensight_dump_test(ut, binary);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstEnsight_Translator.cc
//---------------------------------------------------------------------------//
