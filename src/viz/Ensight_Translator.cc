//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Translator.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator implementation file (non-templated code).
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Ensight_Translator.hh"

namespace rtt_viz
{

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for Ensight_Translator.
 *
 * The constructor automatically knows the number of Ensight cell types (7),
 * and it sets data dependent upon this variable appropriately.
 *
 * \param num_dumps number of timesteps dumped to Ensight 
 */
Ensight_Translator::Ensight_Translator(int num_dumps)
    : num_ensight_cell_types(15),
      ensight_cell_names(num_ensight_cell_types),
      vrtx_cnt(num_ensight_cell_types),
      cell_type_index(num_ensight_cell_types),
      ntime_current(0),
      dump_times(num_dumps, 0.0),
      igrdump_num(0)
{
    // assign values to ensight_cell_names
    ensight_cell_names[0]  = "point";
    ensight_cell_names[1]  = "bar2";
    ensight_cell_names[2]  = "bar3";
    ensight_cell_names[3]  = "tria3";
    ensight_cell_names[4]  = "tria6";
    ensight_cell_names[5]  = "quad4";
    ensight_cell_names[6]  = "quad8";
    ensight_cell_names[7]  = "tetra4";
    ensight_cell_names[8]  = "tetra10";
    ensight_cell_names[9]  = "pyramid5";
    ensight_cell_names[10] = "pyramid13";
    ensight_cell_names[11] = "hexa8";
    ensight_cell_names[12] = "hexa20";
    ensight_cell_names[13] = "penta6";
    ensight_cell_names[14] = "penta15";

    // assign values to vrtx_center
    vrtx_cnt[0]  = 1;
    vrtx_cnt[1]  = 2;
    vrtx_cnt[2]  = 3;
    vrtx_cnt[3]  = 3;
    vrtx_cnt[4]  = 6;
    vrtx_cnt[5]  = 4;
    vrtx_cnt[6]  = 8;
    vrtx_cnt[7]  = 4;
    vrtx_cnt[8]  = 10;
    vrtx_cnt[9]  = 5;
    vrtx_cnt[10] = 13;
    vrtx_cnt[11] = 8;
    vrtx_cnt[12] = 20;
    vrtx_cnt[13] = 6;
    vrtx_cnt[14] = 15;

    // assign values to cell_type_index
    cell_type_index[0]  = point;
    cell_type_index[1]  = two_node_bar;
    cell_type_index[2]  = three_node_bar;
    cell_type_index[3]  = three_node_triangle;
    cell_type_index[4]  = six_node_triangle;
    cell_type_index[5]  = four_node_quadrangle;
    cell_type_index[6]  = eight_node_quadrangle;
    cell_type_index[7]  = four_node_tetrahedron;
    cell_type_index[8]  = ten_node_tetrahedron;
    cell_type_index[9]  = five_node_pyramid;
    cell_type_index[10] = thirteen_node_pyramid;
    cell_type_index[11] = eight_node_hexahedron;
    cell_type_index[12] = twenty_node_hexahedron;
    cell_type_index[13] = six_node_wedge;
    cell_type_index[14] = fifteen_node_wedge;
}

} // end of rtt_viz

//---------------------------------------------------------------------------//
//                              end of Ensight_Translator.cc
//---------------------------------------------------------------------------//
