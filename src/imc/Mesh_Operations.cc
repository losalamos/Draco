//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Mesh_Operations.cc
 * \author Thomas M. Evans
 * \date   Mon Dec 20 15:59:32 1999
 * \brief  Mesh_Operations implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Mesh_Operations.hh"
#include "c4/global.hh"
#include <cmath>

namespace rtt_imc
{

using rtt_mc::OS_Mesh;
using rtt_mc::Coord_sys;
using rtt_rng::Sprng;
using dsxx::SP;
using std::vector;

//===========================================================================//
// OS_MESH SPECIALIZATION IMPLEMENTATION AND INTERFACE
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Mesh_Operations<OS_Mesh>::Mesh_Operations(SP_Mesh mesh, SP_Mat_State state,
					  SP_Topology topology)
    : t4_slope(mesh)
{
    Require(mesh->num_cells() == topology->num_cells(C4::node()));

    // build the T^4 data based upon the topology
    if (topology->get_parallel_scheme() == "replication")
	build_replication_T4_slope(state);
    else if (topology->get_parallel_scheme() == "DD")
	build_DD_T4_slope(state, topology);
    else
    {
	Insist(0, "We don't have general support yet!");
    }
}

//---------------------------------------------------------------------------//
// SAMPLE PARTICLE POSITION WITH A TILT
//---------------------------------------------------------------------------//

Mesh_Operations<OS_Mesh>::sf_double 
Mesh_Operations<OS_Mesh>::sample_pos_tilt(int cell, 
					  double T, 
					  Sprng &random) const
{    
    // set coord system and mesh
    const OS_Mesh &mesh = t4_slope.get_Mesh();
    SP<Coord_sys> coord = t4_slope.get_Mesh().get_SPCoord();

    // assign minimums and maximums for cells dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = mesh.min(d, cell);
	vmax[d-1] = mesh.max(d, cell);
    }

    // use coord_sys to sample the location
    double T4             = std::pow(T, 4);
    vector<double> slopes = t4_slope(cell);
    vector<double> r = 
	coord->sample_pos(vmin, vmax, random, slopes, T4);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// BUILD T4 SLOPES
//---------------------------------------------------------------------------//
// Build the T4 slopes in a full replication topology

void Mesh_Operations<OS_Mesh>::build_replication_T4_slope(SP_Mat_State state)
{
    Require(state->num_cells() == t4_slope.get_Mesh().num_cells());

    // set number of cells
    int num_cells = state->num_cells();

    // get a reference to the mesh
    const OS_Mesh &mesh = t4_slope.get_Mesh();

    // T4 values used in calculation
    double t4_low;
    double t4_high;
    double delta_r;

    // sweep through cells and build T4 slopes
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	// calculate 4th power of cell temperature
	double t4 = std::pow(state->get_T(cell), 4);

	// sweep through coordinates
	for ( int coord = 1; coord <= t4_slope.size(); coord++)
	{
	    Check(t4_slope.size(coord) == num_cells);

	    // get face indices --> these face indices are OS_Mesh dependent
	    int face_low  = 2*coord - 1;
	    int face_high = 2*coord;
	    int cell_low  = mesh.next_cell(cell, face_low);
	    int cell_high = mesh.next_cell(cell, face_high);

	    // set slope to zero if either side is radiatively reflecting
	    if (cell_low == cell || cell_high == cell)
		t4_slope(coord, cell) = 0.0;

	    // set slope to zero if both sides are radiatively vacuum
	    else if (cell_low == 0 && cell_high == 0)
		t4_slope(coord, cell) = 0.0;

	    // if low side is vacuum, use only two t^4's
	    else if (cell_low == 0)
	    {
		t4_high = std::pow(state->get_T(cell_high), 4);
		delta_r = 0.5 * (mesh.dim(coord, cell) + 
				 mesh.dim(coord, cell_high));

		t4_slope(coord, cell) = (t4_high - t4) / delta_r;

		// make sure slope isn't too large so as to give a negative
		// t4_low.  If so, limit slope so t4_low is zero.
		t4_low = t4 - t4_slope(coord, cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_low < 0.0)
		    t4_slope(coord, cell) = 2.0 * t4 / mesh.dim(coord, cell); 
	    }

	    // if high side is vacuum, use only two t^4's
	    else if (cell_high == 0)
	    {
		t4_low = std::pow(state->get_T(cell_low), 4);
		delta_r = 0.5 * (mesh.dim(coord, cell) + 
				 mesh.dim(coord, cell_low));
		t4_slope(coord, cell) = (t4 - t4_low) / delta_r;

		// make sure slope isn't too large so as to give a negative
		// t4_high.  If so, limit slope so t4_high is zero.
		t4_high = t4 + t4_slope(coord,cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = 2.0 * t4 / mesh.dim(coord, cell);
	    }

	    // no conditions on calculating slope; just do it
	    else
	    {
		t4_low  = std::pow(state->get_T(cell_low),  4);
		t4_high = std::pow(state->get_T(cell_high), 4);

		double low_slope = (t4 - t4_low) /
		    (0.5 * (mesh.dim(coord, cell_low) +
			    mesh.dim(coord, cell)) );

		double high_slope = (t4_high - t4) /
		    (0.5 * (mesh.dim(coord, cell) +
			    mesh.dim(coord, cell_high)) );

		double t4_lo_edge = t4 - low_slope  * 0.5 * 
		    mesh.dim(coord, cell);
		double t4_hi_edge = t4 + high_slope * 0.5 *
		    mesh.dim(coord, cell);

		t4_slope(coord, cell) = (t4_hi_edge - t4_lo_edge) / 
		    mesh.dim(coord, cell);
	    }
	}
    }
}

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                              end of Mesh_Operations.cc
//---------------------------------------------------------------------------//
