//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/DD_Mesh.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 26 13:26:21 2000
 * \brief  DD_Mesh header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_test_DD_Mesh_hh__
#define __imc_test_DD_Mesh_hh__

#include "../Mat_State.hh"
#include "mc/Coord_sys.hh"
#include "mc/XYCoord_sys.hh"
#include "mc/Layout.hh"
#include "mc/OS_Mesh.hh"
#include "mc/General_Topology.hh"
#include "c4/global.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>

namespace rtt_imc_test
{

//===========================================================================//
// BUILD DD MESHES
//===========================================================================//
// build 9 cell meshes 2 cells on proc 0; 2 cells on proc 1; 2 cells on proc
// 2; 3 cells on proc 3, +y side is reflecting

rtt_dsxx::SP<rtt_mc::OS_Mesh> build_Mesh() 
{
    using rtt_mc::OS_Mesh;
    using rtt_mc::XYCoord_sys;
    using rtt_mc::Layout;
    using rtt_mc::Coord_sys;
    using rtt_dsxx::SP;

    Require (C4::nodes() == 4);

    SP<OS_Mesh> mesh;

    if (C4::node() == 0)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(2);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = 0;
	lay(1, 2) = 2;
	lay(1, 3) = 0;
	lay(1, 4) = -2;

	lay(2, 1) = 1;
	lay(2, 2) = -1;
	lay(2, 3) = 0;
	lay(2, 4) = -3;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(6);
	vertex[1].resize(6);
	
	vertex[0][0] = 0;
	vertex[0][1] = 1;
	vertex[0][2] = 2;
	vertex[0][3] = 0;
	vertex[0][4] = 1;
	vertex[0][5] = 2;
	vertex[1][0] = 0;
	vertex[1][1] = 0;
	vertex[1][2] = 0;
	vertex[1][3] = 1;
	vertex[1][4] = 1;
	vertex[1][5] = 1;

	OS_Mesh::vf_int cell_pair(2);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 4;
	cell_pair[0][3] = 5;
	cell_pair[1][0] = 2;
	cell_pair[1][1] = 3;
	cell_pair[1][2] = 5;
	cell_pair[1][3] = 6;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }
    else if (C4::node() == 1)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(2);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = -2;
	lay(1, 2) = 0;
	lay(1, 3) = 0;
	lay(1, 4) = -4;

	lay(2, 1) = 0;
	lay(2, 2) = -3;
	lay(2, 3) = -1;
	lay(2, 4) = -5;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(8);
	vertex[1].resize(8);
	
	vertex[0][0] = 2;
	vertex[0][1] = 3;
	vertex[0][2] = 0;
	vertex[0][3] = 1;
	vertex[0][4] = 2;
	vertex[0][5] = 3;
	vertex[0][6] = 0;
	vertex[0][7] = 1;
	vertex[1][0] = 0;
	vertex[1][1] = 0;
	vertex[1][2] = 1;
	vertex[1][3] = 1;
	vertex[1][4] = 1;
	vertex[1][5] = 1;
	vertex[1][6] = 2;
	vertex[1][7] = 2;

	OS_Mesh::vf_int cell_pair(2);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 5;
	cell_pair[0][3] = 6;
	cell_pair[1][0] = 3;
	cell_pair[1][1] = 4;
	cell_pair[1][2] = 7;
	cell_pair[1][3] = 8;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }
    else if (C4::node() == 2)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(2);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = -3;
	lay(1, 2) = 2;
	lay(1, 3) = -1;
	lay(1, 4) = -4;

	lay(2, 1) = 1;
	lay(2, 2) = 0;
	lay(2, 3) = -2;
	lay(2, 4) = -5;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(6);
	vertex[1].resize(6);
	
	vertex[0][0] = 1;
	vertex[0][1] = 2;
	vertex[0][2] = 3;
	vertex[0][3] = 1;
	vertex[0][4] = 2;
	vertex[0][5] = 3;
	vertex[1][0] = 1;
	vertex[1][1] = 1;
	vertex[1][2] = 1;
	vertex[1][3] = 2;
	vertex[1][4] = 2;
	vertex[1][5] = 2;

	OS_Mesh::vf_int cell_pair(2);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 4;
	cell_pair[0][3] = 5;
	cell_pair[1][0] = 2;
	cell_pair[1][1] = 3;
	cell_pair[1][2] = 5;
	cell_pair[1][3] = 6;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }
    else if (C4::node() == 3)
    {
	SP<Coord_sys> cs(new XYCoord_sys());

	// build layout
	Layout lay(3);
	for (int i = 1; i <= lay.num_cells(); i++)
	    lay.set_size(i, 4);

	lay(1, 1) = 0;
	lay(1, 2) = 2;
	lay(1, 3) = -1;
	lay(1, 4) = 1;

	lay(2, 1) = 1;
	lay(2, 2) = 3;
	lay(2, 3) = -2;
	lay(2, 4) = 2;

	lay(3, 1) = 2;
	lay(3, 2) = 0;
	lay(3, 3) = -3;
	lay(3, 4) = 3;

	// build coordinate data
	OS_Mesh::vf_double vertex(2);
	vertex[0].resize(8);
	vertex[1].resize(8);
	
	vertex[0][0] = 0;
	vertex[0][1] = 1;
	vertex[0][2] = 2;
	vertex[0][3] = 3;
	vertex[0][4] = 0;
	vertex[0][5] = 1;
	vertex[0][6] = 2;
	vertex[0][7] = 3;
	vertex[1][0] = 2;
	vertex[1][1] = 2;
	vertex[1][2] = 2;
	vertex[1][3] = 2;
	vertex[1][4] = 3;
	vertex[1][5] = 3;
	vertex[1][6] = 3;
	vertex[1][7] = 3;

	OS_Mesh::vf_int cell_pair(3);
	cell_pair[0].resize(4);
	cell_pair[1].resize(4);
	cell_pair[2].resize(4);
	
	cell_pair[0][0] = 1;
	cell_pair[0][1] = 2;
	cell_pair[0][2] = 5;
	cell_pair[0][3] = 6;
	cell_pair[1][0] = 2;
	cell_pair[1][1] = 3;
	cell_pair[1][2] = 6;
	cell_pair[1][3] = 7;
	cell_pair[2][0] = 3;
	cell_pair[2][1] = 4;
	cell_pair[2][2] = 7;
	cell_pair[2][3] = 8;
	
	mesh = new OS_Mesh(cs, lay, vertex, cell_pair, true);
    }

    return mesh;
}

//===========================================================================//
// BUILD DD TOPOLOGIES
//===========================================================================//
// build a General Topology for 9 cell mesh described above

rtt_dsxx::SP<rtt_mc::Topology> build_Topology()
{
    using rtt_mc::Topology;
    using rtt_mc::General_Topology;
    using rtt_dsxx::SP;
    using std::vector;

    Require (C4::nodes() == 4);

    Topology::vf_int cpp(4);
    Topology::vf_int ppc(9, vector<int>(1));
    Topology::vf_int bc(4);

    // cells-per-proc data
    cpp[0].resize(2);
    cpp[0][0] = 1;
    cpp[0][1] = 2;

    cpp[1].resize(2);
    cpp[1][0] = 3;
    cpp[1][1] = 4;

    cpp[2].resize(2);
    cpp[2][0] = 5;
    cpp[2][1] = 6;

    cpp[3].resize(3);
    cpp[3][0] = 7;
    cpp[3][1] = 8;
    cpp[3][2] = 9;

    // procs_per_cell data
    ppc[0][0] = 0;
    ppc[1][0] = 0;
    ppc[2][0] = 1;
    ppc[3][0] = 1;
    ppc[4][0] = 2;
    ppc[5][0] = 2;
    ppc[6][0] = 3;
    ppc[7][0] = 3;
    ppc[8][0] = 3;

    // boundary cell data
    bc[0].resize(3);
    bc[0][0] = 3;
    bc[0][1] = 4;
    bc[0][2] = 5;

    bc[1].resize(5);
    bc[1][0] = 1;
    bc[1][1] = 2;
    bc[1][2] = 5;
    bc[1][3] = 6;
    bc[1][4] = 7;

    bc[2].resize(5);
    bc[2][0] = 2;
    bc[2][1] = 3;
    bc[2][2] = 4;
    bc[2][3] = 8;
    bc[2][4] = 9;
    
    bc[3].resize(3);
    bc[3][0] = 4;
    bc[3][1] = 5;
    bc[3][2] = 6;

    SP<Topology> topology(new General_Topology(cpp, ppc, bc, "DD"));
    return topology;
}

//===========================================================================//
// BUILD DD MAT_STATES
//===========================================================================//
// build mat states with temperatures = C4::node() + 10.5

rtt_dsxx::SP<rtt_imc::Mat_State<rtt_mc::OS_Mesh> >
build_Mat(rtt_dsxx::SP<rtt_mc::OS_Mesh> mesh)
{
    using rtt_imc::Mat_State;
    using rtt_mc::OS_Mesh;
    using rtt_dsxx::SP;

    Require (C4::nodes() == 4);

    OS_Mesh::CCSF_double x(mesh);
    OS_Mesh::CCSF_double temp(mesh);

    for (int i = 1; i <= temp.size(); i++)
	temp(i) = C4::node() + 10.5;
    
    SP<Mat_State<OS_Mesh> > mat(new Mat_State<OS_Mesh>(x,temp,x,x));
    return mat;
}

} // end namespace rtt_imc_test

#endif                          // __imc_test_DD_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/DD_Mesh.hh
//---------------------------------------------------------------------------//
