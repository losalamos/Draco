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
#include "mc/Parallel_Data_Operator.hh"
#include "mc/AMR_Layout.hh"
#include "c4/global.hh"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace rtt_imc
{

//===========================================================================//
// OS_MESH SPECIALIZATION IMPLEMENTATION AND INTERFACE
//===========================================================================//

using rtt_mc::OS_Mesh;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Mesh_Operations<OS_Mesh>::Mesh_Operations(SP_Mesh mesh, 
					  SP_Mat_State state,
					  SP_Topology topology,
					  SP_Comm_Patterns comm_patterns)
    : t4_slope(mesh)
{
    Require (mesh);
    Require (state);
    Require (topology);
    Require (comm_patterns);
    Require (mesh->num_cells() == topology->num_cells(C4::node()));

    // build the T^4 data based upon the topology
    if (topology->get_parallel_scheme() == "replication")
    {
	Check (!(*comm_patterns));
	build_replication_T4_slope(state);
    }
    else if (topology->get_parallel_scheme() == "DD")
    {
	Check (*comm_patterns);
	build_DD_T4_slope(state, topology, comm_patterns);
    }
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
					  rtt_rng::Sprng &random) const
{    
    using rtt_mc::Coord_sys;
    using rtt_dsxx::SP;
    using std::vector;

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
    double T4            = std::pow(T, 4);
    vector<double> slope = t4_slope(cell);
    vector<double> r     = coord->sample_pos(vmin, vmax, random, slope, T4);

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
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell);
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

		// put checks to make sure that slope on high and low ends
		// are not too large
		t4_high = t4 + t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		t4_low  = t4 - t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		Check (t4_high >= 0.0 || t4_low >= 0.0);

		// if high is less than 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell); 
		else if (t4_low < 0.0)
		    t4_slope(coord, cell) = 2. * t4 / mesh.dim(coord, cell);
	    }
	}
    }
}

//---------------------------------------------------------------------------//
// build the T4 slopes in a full DD topology

void Mesh_Operations<OS_Mesh>::build_DD_T4_slope(SP_Mat_State state,
						 SP_Topology topology,
						 SP_Comm_Patterns com_pat)
{
    Require(state->num_cells() == t4_slope.get_Mesh().num_cells());

    // define the boundary cell fields that we need to perform the
    // calculation 
    sf_double bc_temp;
    vf_double bc_dim;;
    
    int num_cells = topology->num_cells(C4::node());
    Check (num_cells == state->num_cells());

    // reference to the mesh
    const OS_Mesh &mesh = t4_slope.get_Mesh();

    // make a Parallel_Data_Operator for performing gathers on boundary cells 
    // with local data arrays
    rtt_mc::Parallel_Data_Operator par_op(topology);

    // get the temperatures on each boundary cell
    {
	// local temperatures
	sf_double local_temp(num_cells);
	for (int cell = 1; cell <= num_cells; cell++)
	    local_temp[cell-1] = state->get_T(cell);

	// fill the boundary cells with temperatures
	par_op.gather_bnd_cell_data(com_pat, local_temp, bc_temp);
    }
    Check (bc_temp.size() == topology->get_boundary_cells(C4::node()));

    // get the widths of each boundary cell
    {
	// size the bc_dim vector to the number of dimensions in this problem
	bc_dim.resize(t4_slope.size());

	// loop through coordinates and fill up width data
	for (int coord = 1; coord <= bc_dim.size(); coord++)
	{
	    // local widths
	    sf_double local_dim(num_cells);
	    for (int cell = 1; cell <= num_cells; cell++)
		local_dim[cell-1] = mesh.dim(coord, cell);

	    // fill up the boundary cells with dimension data
	    par_op.gather_bnd_cell_data(com_pat, local_dim, bc_dim[coord-1]);
	    Check (bc_dim[coord-1].size() ==
		   topology->get_boundary_cells(C4::node()));
	}
    }
    Check (bc_dim.size() == mesh.get_Coord().get_dim());

    // T4 values used in calculation
    double t4_low;
    double t4_high;
    double dim_low;
    double dim_high;
    double delta_r;

    // sweep through cells and build T4 slopes
    for (int cell = 1; cell <= num_cells; cell++)
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

	    // calculate "vanilla" t4 high and low side values, in
	    // particular, if the cell has a negative number then it is a
	    // boundary cell and we need to get its temperature value from
	    // the boundary temperature field
	    if (cell_low < 0)
	    {
		t4_low  = std::pow(bc_temp[-cell_low - 1], 4);
		dim_low = bc_dim[coord-1][-cell_low - 1];
	    }
	    else if (cell_low > 0)
	    {
		t4_low  = std::pow(state->get_T(cell_low), 4);
		dim_low = mesh.dim(coord, cell_low); 
	    }

	    if (cell_high < 0)
	    {
		t4_high  = std::pow(bc_temp[-cell_high - 1], 4);
		dim_high = bc_dim[coord-1][-cell_high - 1];
	    }
	    else if (cell_high > 0)
	    {
		t4_high  = std::pow(state->get_T(cell_high), 4);
		dim_high = mesh.dim(coord, cell_high);
	    }

	    // calculate the T4 slopes

	    // set slope to zero if either side is radiatively reflecting
	    if (cell_low == cell || cell_high == cell)
		t4_slope(coord, cell) = 0.0;

	    // set slope to zero if both sides are radiatively vacuum
	    else if (cell_low == 0 && cell_high == 0)
		t4_slope(coord, cell) = 0.0;

	    // if low side is vacuum, use only two t^4's
	    else if (cell_low == 0)
	    {
		delta_r = 0.5 * (mesh.dim(coord, cell) + dim_high);

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
		delta_r = 0.5 * (mesh.dim(coord, cell) + dim_low);
		t4_slope(coord, cell) = (t4 - t4_low) / delta_r;

		// make sure slope isn't too large so as to give a negative
		// t4_high.  If so, limit slope so t4_high is zero.
		t4_high = t4 + t4_slope(coord,cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell);
	    }

	    // no conditions on calculating slope; just do it
	    else
	    {
		double low_slope = (t4 - t4_low) /
		    (0.5 * (dim_low + mesh.dim(coord, cell)) );

		double high_slope = (t4_high - t4) /
		    (0.5 * (mesh.dim(coord, cell) + dim_high) );

		double t4_lo_edge = t4 - low_slope  * 0.5 * 
		    mesh.dim(coord, cell); 
		double t4_hi_edge = t4 + high_slope * 0.5 *
		    mesh.dim(coord, cell);

		t4_slope(coord, cell) = (t4_hi_edge - t4_lo_edge) / 
		    mesh.dim(coord, cell);

		// put checks to make sure that slope on high and low ends
		// are not too large
		t4_high = t4 + t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		t4_low  = t4 - t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		Check (t4_high >= 0.0 || t4_low >= 0.0);

		// if high is less than 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell); 
		else if (t4_low < 0.0)
		    t4_slope(coord, cell) = 2. * t4 / mesh.dim(coord, cell);
	    }
	}
    }
}

//===========================================================================//
// RZWEDGE_MESH SPECIALIZATIONS
//===========================================================================//

using rtt_mc::RZWedge_Mesh;
using rtt_mc::AMR_Layout;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Mesh_Operations<RZWedge_Mesh>::Mesh_Operations(SP_Mesh mesh,
					       SP_Mat_State state,
					       SP_Topology topology,
					       SP_Comm_Patterns patterns)
    : t4_slope(mesh)
{
    Require (mesh);
    Require (state);
    Require (topology);
    Require (patterns);
    Require (mesh->num_cells() == topology->num_cells(C4::node()));

    // RZWedge_Meshes always have 3 dimensions
    Check(t4_slope.size() == 3);

    // build the T^4 data based upon the topology
    if (topology->get_parallel_scheme() == "replication")
    {
	Check (!(*patterns));
	build_replication_T4_slope(state);
    }
    else if (topology->get_parallel_scheme() == "DD")
    {
	Check (*patterns);
	build_DD_T4_slope(state, topology, patterns);
    }
    else
    {
	Insist(0, "We don't have general support yet!");
    }
}

//---------------------------------------------------------------------------//
// SAMPLE PARTICLE POSITION WITH A TILT
//---------------------------------------------------------------------------//

Mesh_Operations<RZWedge_Mesh>::sf_double 
Mesh_Operations<RZWedge_Mesh>::sample_pos_tilt(int cell, 
					       double T, 
					       rtt_rng::Sprng &random) const
{    
    using std::vector;

    // set coord system and mesh
    const RZWedge_Mesh &mesh = t4_slope.get_Mesh();
    Check (mesh.get_SPCoord()->get_Coord() == "xyz");

    // return position
    vector<double> r;

    // T4 slopes and T4 temperature
    double T4            = std::pow(T, 4);
    vector<double> slope = t4_slope(cell);

    // for now sample uniformly
    r = mesh.sample_pos(cell, random, slope, T4);
    Check (r.size() == 3);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// BUILD T4 SLOPES
//---------------------------------------------------------------------------//
// Build the T4 slopes in a full replication topology

void Mesh_Operations<RZWedge_Mesh>::build_replication_T4_slope(SP_Mat_State 
							       mat_state)
{
    Require(mat_state->num_cells() == t4_slope.get_Mesh().num_cells());

    // set number of cells
    int num_cells = mat_state->num_cells();

    // get a reference to the mesh and layout
    const RZWedge_Mesh &mesh = t4_slope.get_Mesh();
    const AMR_Layout &layout = mesh.get_Layout();

    // T4 values used in calculation
    double t4_low  = 0;
    double t4_high = 0;
    double delta_r = 0;

    // sweep through cells and build T4 slopes
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	// calculate 4th power of cell temperature
	double t4 = std::pow(mat_state->get_T(cell), 4);

	// sweep through coordinates
	for ( int coord = 1; coord <= 3; coord += 2 /* skip y */)
	{
	    Check(t4_slope.size(coord) == num_cells);

	    // get coarse face indices
	    int coarse_face_low  = 2*coord - 1;
	    int coarse_face_high = 2*coord;

	    // determine number of cells across a face
	    int num_across_low  = layout.num_cells_across(cell,
							  coarse_face_low);
	    int num_across_high = layout.num_cells_across(cell,
							  coarse_face_high);

	    // get t4_high and t4_low
	    int cell_low;
	    int cell_high;
	    
	    // T4 low
	    t4_low = 0.0;
	    for (int i = 1; i <= num_across_low; i++)
	    {
		cell_low = layout(cell, coarse_face_low, i);
		if (cell_low > 0)
		    t4_low  += std::pow(mat_state->get_T(cell_low), 4);
	    }
	    t4_low = t4_low / static_cast<double>(num_across_low);

	    // T4 high
	    t4_high = 0.0;
	    for (int i = 1; i <= num_across_high; i++)
	    {
		cell_high = layout(cell, coarse_face_high, i);
		if (cell_high > 0)
		    t4_high  += std::pow(mat_state->get_T(cell_high), 4);
	    }
	    t4_high = t4_high / static_cast<double>(num_across_high);

	    // explicitly set cell lows and cell highs
	    cell_low  = layout(cell, coarse_face_low, 1);
	    cell_high = layout(cell, coarse_face_high, 1);

	    // set slope to zero if either side is radiatively reflecting
	    if (cell_low == cell || cell_high == cell)
		t4_slope(coord, cell) = 0.0;

	    // set slope to zero if both sides are radiatively vacuum
	    else if (cell_low == 0 && cell_high == 0)
		t4_slope(coord, cell) = 0.0;

	    // if low side is vacuum, use only two t^4's
	    else if (cell_low == 0)
	    {
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
		delta_r = 0.5 * (mesh.dim(coord, cell) + 
				 mesh.dim(coord, cell_low));
		t4_slope(coord, cell) = (t4 - t4_low) / delta_r;

		// make sure slope isn't too large so as to give a negative
		// t4_high.  If so, limit slope so t4_high is zero.
		t4_high = t4 + t4_slope(coord,cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell);
	    }

	    // no conditions on calculating slope; just do it
	    else
	    {
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

		// put checks to make sure that slope on high and low ends
		// are not too large
		t4_high = t4 + t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		t4_low  = t4 - t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		Check (t4_high >= 0.0 || t4_low >= 0.0);

		// if high is less than 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell); 
		else if (t4_low < 0.0)
		    t4_slope(coord, cell) = 2. * t4 / mesh.dim(coord, cell);
	    }
	}

	// explicitly set y coordinate t4_slopes to zero
	t4_slope(2, cell) = 0.0;
    }
}

//---------------------------------------------------------------------------//
// build the T4 slopes in a full DD topology

void Mesh_Operations<RZWedge_Mesh>::build_DD_T4_slope(SP_Mat_State mat_state,
						      SP_Topology topology,
						      SP_Comm_Patterns com_pat)
{
    Require(mat_state->num_cells() == t4_slope.get_Mesh().num_cells());

    // define the boundary cell fields that we need to perform the
    // calculation
    sf_double bc_temp;
    vf_double bc_dim;

    // define the number of local cells
    int num_cells = topology->num_cells(C4::node());
    Check (num_cells == mat_state->num_cells());

    // get a reference to the mesh and layout
    const RZWedge_Mesh &mesh = t4_slope.get_Mesh();
    const AMR_Layout &layout = mesh.get_Layout();

    // make a Parallel_Data_Operator for performing gathers on boundary cells 
    // with local data arrays
    rtt_mc::Parallel_Data_Operator par_op(topology);

    // get the temperatures on each boundary cell
    {
	// local temperatures
	sf_double local_temp(num_cells);
	for (int cell = 1; cell <= num_cells; cell++)
	    local_temp[cell-1] = mat_state->get_T(cell);

	// fill the boundary cells with temperatures
	par_op.gather_bnd_cell_data(com_pat, local_temp, bc_temp);
    }
    Check (bc_temp.size() == topology->get_boundary_cells(C4::node()));

    // get the widths of each boundary cell
    {
	// size the bc_dim vector to the number of dimensions in this problem 
	bc_dim.resize(3);

	// loop through coordinates and fill up width data
	for (int coord = 1; coord <= 3; coord += 2)
	{
	    // local widths
	    sf_double local_dim(num_cells);
	    for (int cell = 1; cell <= num_cells; cell++)
		local_dim[cell-1] = mesh.dim(coord, cell);

	    // fill up the boundary cells with dimension data
	    par_op.gather_bnd_cell_data(com_pat, local_dim, bc_dim[coord-1]);
	    Check (bc_dim[coord-1].size() ==
		   topology->get_boundary_cells(C4::node()));
	}

	// fill y dimensions with 0 width
	bc_dim[1].resize(topology->get_boundary_cells(C4::node()));
	std::fill(bc_dim[1].begin(), bc_dim[1].end(), 0.0);
    }
    Check (bc_dim.size() == mesh.get_Coord().get_dim());

    // T4 values used in calculation
    double t4_low   = 0;
    double t4_high  = 0;
    double dim_low  = 0;
    double dim_high = 0;
    double delta_r  = 0;

    // sweep through cells and build T4 slopes
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	// calculate 4th power of cell temperature
	double t4 = std::pow(mat_state->get_T(cell), 4);

	// sweep through coordinates
	for ( int coord = 1; coord <= 3; coord += 2 /* skip y */)
	{
	    Check(t4_slope.size(coord) == num_cells);

	    // get coarse face indices
	    int coarse_face_low  = 2*coord - 1;
	    int coarse_face_high = 2*coord;

	    // determine number of cells across a face
	    int num_across_low  = layout.num_cells_across(cell,
							  coarse_face_low);
	    int num_across_high = layout.num_cells_across(cell,
							  coarse_face_high);

	    // get t4_high and t4_low; dim_high and dim_low
	    int cell_low;
	    int cell_high;
	    
	    // T4 low
	    t4_low = 0.0;
	    for (int i = 1; i <= num_across_low; i++)
	    {
		cell_low = layout(cell, coarse_face_low, i);
		if (cell_low > 0)
		{
		    t4_low += std::pow(mat_state->get_T(cell_low), 4);
		    dim_low = mesh.dim(coord, cell_low);
		}
		else if (cell_low < 0)
		{
		    t4_low += std::pow(bc_temp[-cell_low - 1], 4);
		    dim_low = bc_dim[coord-1][-cell_low - 1];
		}
	    }
	    t4_low = t4_low / static_cast<double>(num_across_low);

	    // T4 high
	    t4_high = 0.0;
	    for (int i = 1; i <= num_across_high; i++)
	    {
		cell_high = layout(cell, coarse_face_high, i);
		if (cell_high > 0)
		{
		    t4_high += std::pow(mat_state->get_T(cell_high), 4);
		    dim_high = mesh.dim(coord, cell_high);
		}
		else if (cell_high < 0)
		{
		    t4_high += std::pow(bc_temp[-cell_high - 1], 4);
		    dim_high = bc_dim[coord-1][-cell_high - 1];
		}
	    }
	    t4_high = t4_high / static_cast<double>(num_across_high);

	    // explicitly set cell lows and cell highs
	    cell_low  = layout(cell, coarse_face_low, 1);
	    cell_high = layout(cell, coarse_face_high, 1);

	    // set slope to zero if either side is radiatively reflecting
	    if (cell_low == cell || cell_high == cell)
		t4_slope(coord, cell) = 0.0;

	    // set slope to zero if both sides are radiatively vacuum
	    else if (cell_low == 0 && cell_high == 0)
		t4_slope(coord, cell) = 0.0;

	    // if low side is vacuum, use only two t^4's
	    else if (cell_low == 0)
	    {
		delta_r = 0.5 * (mesh.dim(coord, cell) + dim_high);

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
		delta_r = 0.5 * (mesh.dim(coord, cell) + dim_low);
		t4_slope(coord, cell) = (t4 - t4_low) / delta_r;

		// make sure slope isn't too large so as to give a negative
		// t4_high.  If so, limit slope so t4_high is zero.
		t4_high = t4 + t4_slope(coord,cell) * 0.5 * 
		    mesh.dim(coord, cell); 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell);
	    }

	    // no conditions on calculating slope; just do it
	    else
	    {
		double low_slope = (t4 - t4_low) /
		    (0.5 * (dim_low + mesh.dim(coord, cell)) );

		double high_slope = (t4_high - t4) /
		    (0.5 * (mesh.dim(coord, cell) + dim_high) );

		double t4_lo_edge = t4 - low_slope  * 0.5 * 
		    mesh.dim(coord, cell);
		double t4_hi_edge = t4 + high_slope * 0.5 *
		    mesh.dim(coord, cell);

		t4_slope(coord, cell) = (t4_hi_edge - t4_lo_edge) / 
		    mesh.dim(coord, cell);

		// put checks to make sure that slope on high and low ends
		// are not too large
		t4_high = t4 + t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		t4_low  = t4 - t4_slope(coord, cell) * .5 * 
		    mesh.dim(coord, cell);

		Check (t4_high >= 0.0 || t4_low >= 0.0);

		// if high is less than 
		if (t4_high < 0.0)
		    t4_slope(coord, cell) = -2. * t4 / mesh.dim(coord, cell); 
		else if (t4_low < 0.0)
		    t4_slope(coord, cell) = 2. * t4 / mesh.dim(coord, cell);
	    }
	}

	// explicitly set y coordinate t4_slopes to zero
	t4_slope(2, cell) = 0.0;
    }
}

//===========================================================================//
// SPHYRAMID_MESH SPECIALIZATIONS
//===========================================================================//

using rtt_mc::Sphyramid_Mesh;
using rtt_mc::Layout;

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//

Mesh_Operations<Sphyramid_Mesh>::Mesh_Operations(SP_Mesh mesh,
						 SP_Mat_State state,
						 SP_Topology topology,
						 SP_Comm_Patterns patterns)
    :t4_slope(mesh)
{
    Require (mesh);
    Require (state);
    Require (topology);
    Require (patterns);
    Require (mesh->num_cells() == topology->num_cells(C4::node()));

    // Sphyramid_Mesh always has 3 dimensions
    Check(this->t4_slope.size() == 3);

    // build the T^4 data based upon the topology
    if (topology->get_parallel_scheme() == "replication")
    {
	Check (!(*patterns));
	build_replication_T4_slope(state);
    }
    else if (topology->get_parallel_scheme() == "DD")
    {

	Check (*patterns);
	build_DD_T4_slope(state, topology, patterns);
    }
    else
    {
	Insist(0, "We don't have general support yet!");
    }

    return;
}

//---------------------------------------------------------------------------//
// SAMPLE PARTICLE POSITION WITH A TILT
//---------------------------------------------------------------------------//

Mesh_Operations<Sphyramid_Mesh>::sf_double
Mesh_Operations<Sphyramid_Mesh>::sample_pos_tilt(int cell, 
						 double T, 
						 rtt_rng::Sprng &random) const
{
    using std::vector;
    using std::pow;
    
    // set coord system and mesh
    const Sphyramid_Mesh &mesh = this->t4_slope.get_Mesh();
    Check (mesh.get_SPCoord()->get_Coord() == "xyz");

    // return position
    vector<double> r;
   
    // T4 slopes and T4 temperature
    double T4            = pow(T,4);
    vector<double> slope = this->t4_slope(cell);

    // sample position
    r = mesh.sample_pos(cell, random, slope, T4);
    Check (r.size() == 3);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// BUILD T4 SLOPES
//---------------------------------------------------------------------------//
// Build the T4 slopes in a full replication topology

void Mesh_Operations<Sphyramid_Mesh>::build_replication_T4_slope(SP_Mat_State 
								 mat_state)
{	    
    using std::pow;

    Require(mat_state->num_cells() == this->t4_slope.get_Mesh().num_cells());

    // set number of cells
    int num_cells = mat_state->num_cells();

    // get a reference to the mesh and layout
    const Sphyramid_Mesh &mesh = this->t4_slope.get_Mesh();
    const Layout &layout       = mesh.get_Layout();

    // sweep through cells and build T4 slopes
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	// calculate this cell's value of t4
	double t4 = pow(mat_state->get_T(cell), 4);

	// >>>calculate x direction slopes<<<
	
	// check size of vector
	Check (this->t4_slope.size(1) == num_cells);

	// get face indices
	int face_low  = 1;
	int face_high = 2;

	// get neighboring cell indices
	int cell_low  = layout(cell, face_low);
	int cell_high = layout(cell, face_high);

	// get neighboring cell temperatures
	double t4_low  = pow(mat_state->get_T(cell_low), 4);
	double t4_high = pow(mat_state->get_T(cell_high), 4);

	// explicity set slope to zero if either face is reflective
	if (cell_low == cell || cell_high == cell)
	{
	    this->t4_slope(1,cell) = 0.0;
	}
	// set slope to zero if both sides are radiatively vacuum
	else if (cell_low == 0 && cell_high == 0)
	{
	    this->t4_slope(1,cell) = 0.0;
	}
	// if low side is vacuum, use only two T^4's
	else if (cell_low == 0)
	{
	    double delta = mesh.high_half_width(cell)
		+mesh.low_half_width(cell_high);
	    
	    this->t4_slope(1,cell) = (t4_high-t4)/delta;

	    // make sure slope isn't too large so as to give a negative
	    // t4_low.  If so, limit slope so t4_low is zero
	    t4_low = t4-this->t4_slope(1,cell)*mesh.low_half_width(cell);
	    if (t4_low < 0.0)
	    {
		this->t4_slope(1,cell) = t4/mesh.low_half_width(cell);
	    }
	}
	// if high side is vacuum, use only two T^4's
	else if (cell_high == 0)
	{
	    double delta = mesh.low_half_width(cell)
		+mesh.high_half_width(cell_low);
	    
	    this->t4_slope(1,cell) = (t4-t4_low)/delta;

	    // make sure slope isn't too large so as to give a negative
	    // t4_high.  If so, limit slope so t4_high is zero
	    t4_high = t4+this->t4_slope(1,cell)*mesh.high_half_width(cell);
	    if (t4_high < 0.0)
	    {
		this->t4_slope(1,cell) = -t4/mesh.high_half_width(cell);
	    }
	}    
	// no conditions on calculating slope; just do it
	else
	{
	    double low_slope = (t4-t4_low)/
		(mesh.low_half_width(cell)+mesh.high_half_width(cell_low));

	    double high_slope = (t4_high-t4)/
		(mesh.high_half_width(cell)+mesh.low_half_width(cell_high));

	    double t4_lo_edge = t4-low_slope*mesh.low_half_width(cell);
	    double t4_hi_edge = t4+high_slope*mesh.high_half_width(cell);
	    
	    this->t4_slope(1,cell) = (t4_hi_edge-t4_lo_edge)/
		(mesh.low_half_width(cell)+mesh.high_half_width(cell));
	   
	    // put checks to make sure that slope is not too large
	    t4_high = t4+this->t4_slope(1,cell)*mesh.high_half_width(cell);
	    t4_low  = t4-this->t4_slope(1,cell)*mesh.low_half_width(cell);
	    
	    // at most one edge can be negative
	    Check (t4_high >= 0.0 || t4_low >= 0.0);
	    
	    // if high edge is negative
	    if (t4_high < 0.0)
	    {
		this->t4_slope(1,cell) = -t4/mesh.high_half_width(cell);
	    }
	    else if (t4_low < 0.0)
	    {
		this->t4_slope(1,cell) = t4/mesh.low_half_width(cell);
	    } 
	    
	}	    
	// explicity set y and z slopes to zero
	Check (this->t4_slope.size(2) == num_cells);
	Check (this->t4_slope.size(3) == num_cells);
	this->t4_slope(2, cell) = 0.0;
	this->t4_slope(3, cell) = 0.0;
    }

    return;
}
//---------------------------------------------------------------------------//
// Build the T4 slopes in a full DD topology

void Mesh_Operations<Sphyramid_Mesh>::build_DD_T4_slope(SP_Mat_State mat_state,
							SP_Topology topology,
							SP_Comm_Patterns com_pat)
{	    
    using std::pow;

    Require(mat_state->num_cells() == this->t4_slope.get_Mesh().num_cells());

    //define the boundary cell fields that we need to perform the calculation
    sf_double bc_temp;
    sf_double bc_high_half_width;
    sf_double bc_low_half_width;
    
    // set number of local cells
    int num_cells = topology->num_cells(C4::node());
    Check (num_cells == mat_state->num_cells());

    // get a reference to the mesh and layout
    const Sphyramid_Mesh &mesh = this->t4_slope.get_Mesh();
    const Layout &layout       = mesh.get_Layout();

    // make a Parallel_Data_Operator for performing gathers on boundary cells
    // with local data arrays
    rtt_mc::Parallel_Data_Operator par_op(topology);

    // get the temperatures on each boundary cell
    {
	// local temperatures
	sf_double local_temp(num_cells);
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    local_temp[cell-1] = mat_state->get_T(cell);
	}
	    // fill the boundary cells with temperatures
	    par_op.gather_bnd_cell_data(com_pat, local_temp, bc_temp);
    }
    Check (bc_temp.size() == topology->get_boundary_cells(C4::node()));

    // get the half-widths of each boundary cell
    {
	// local half-widths
	sf_double local_low(num_cells);
	sf_double local_high(num_cells);
	for (int cell = 1; cell <= num_cells; cell++)
	{
	    local_low[cell-1]  = mesh.low_half_width(cell);
	    local_high[cell-1] = mesh.high_half_width(cell);
	}

	// fill up the boundary cells with width data
	par_op.gather_bnd_cell_data(com_pat, local_low, bc_low_half_width);
	par_op.gather_bnd_cell_data(com_pat, local_high, bc_high_half_width);

    }

    // sweep through cells and build T4 slopes
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	// calculate this cell's value of t4
	double t4 = pow(mat_state->get_T(cell), 4);

	// >>>calculate x direction slopes<<<
	
	// check size of vector
	Check (this->t4_slope.size(1) == num_cells);

	// get face indices
	int face_low  = 1;
	int face_high = 2;

	// get neighboring cell indices
	int cell_low  = layout(cell, face_low);
	int cell_high = layout(cell, face_high);

	// get neighboring cell temperatures and half_widhts
	double t4_low                = 0.0;
	double t4_high               = 0.0;
	double lo_cell_hi_half_width = 0.0;
	double hi_cell_lo_half_width = 0.0;
	if (cell_low > 0)
	{
	    t4_low                = pow(mat_state->get_T(cell_low), 4);
	    lo_cell_hi_half_width = mesh.high_half_width(cell_low);
	}
	else if (cell_low < 0)
	{
	    t4_low                = pow(bc_temp[-cell_low-1], 4);
	    lo_cell_hi_half_width = bc_high_half_width[-cell_low-1];
	}
	if (cell_high > 0)
	{
	    t4_high               = pow(mat_state->get_T(cell_high), 4);
	    hi_cell_lo_half_width = mesh.low_half_width(cell_high);
	}
	else if (cell_high < 0)
	{
	    t4_high               = pow(bc_temp[-cell_high-1], 4);
	    hi_cell_lo_half_width = bc_low_half_width[-cell_high-1];
	}

	// explicity set slope to zero if either face is reflective
	if (cell_low == cell || cell_high == cell)
	{
	    this->t4_slope(1,cell) = 0.0;
	}
	// set slope to zero if both sides are radiatively vacuum
	else if (cell_low == 0 && cell_high == 0)
	{
	    this->t4_slope(1,cell) = 0.0;
	}
	// if low side is vacuum, use only two T^4's
	else if (cell_low == 0)
	{
	    double delta = mesh.high_half_width(cell)
		+hi_cell_lo_half_width;
	    
	    this->t4_slope(1,cell) = (t4_high-t4)/delta;

	    // make sure slope isn't too large so as to give a negative
	    // t4_low.  If so, limit slope so t4_low is zero
	    t4_low = t4-this->t4_slope(1,cell)*mesh.low_half_width(cell);
	    if (t4_low < 0.0)
	    {
		this->t4_slope(1,cell) = t4/mesh.low_half_width(cell);
	    }
	}
	// if high side is vacuum, use only two T^4's
	else if (cell_high == 0)
	{
	    double delta = mesh.low_half_width(cell)
		+lo_cell_hi_half_width;
	    
	    this->t4_slope(1,cell) = (t4-t4_low)/delta;

	    // make sure slope isn't too large so as to give a negative
	    // t4_high.  If so, limit slope so t4_high is zero
	    t4_high = t4+this->t4_slope(1,cell)*mesh.high_half_width(cell);
	    if (t4_high < 0.0)
	    {
		this->t4_slope(1,cell) = -t4/mesh.high_half_width(cell);
	    }
	}    
	// no conditions on calculating slope; just do it
	else
	{
	    double low_slope = (t4-t4_low)/
		(mesh.low_half_width(cell)+lo_cell_hi_half_width);

	    double high_slope = (t4_high-t4)/
		(mesh.high_half_width(cell)+hi_cell_lo_half_width);

	    double t4_lo_edge = t4-low_slope*mesh.low_half_width(cell);
	    double t4_hi_edge = t4+high_slope*mesh.high_half_width(cell);
	    
	    this->t4_slope(1,cell) = (t4_hi_edge-t4_lo_edge)/
		(mesh.low_half_width(cell)+mesh.high_half_width(cell));
	   
	    // put checks to make sure that slope is not too large
	    t4_high = t4+this->t4_slope(1,cell)*mesh.high_half_width(cell);
	    t4_low  = t4-this->t4_slope(1,cell)*mesh.low_half_width(cell);
	    
	    // at most one edge can be negative
	    Check (t4_high >= 0.0 || t4_low >= 0.0);
	    
	    // if high edge is negative
	    if (t4_high < 0.0)
	    {
		this->t4_slope(1,cell) = -t4/mesh.high_half_width(cell);
	    }
	    else if (t4_low < 0.0)
	    {
		this->t4_slope(1,cell) = t4/mesh.low_half_width(cell);
	    } 
	    
	}	    
	// explicity set y and z slopes to zero
	Check (this->t4_slope.size(2) == num_cells);
	Check (this->t4_slope.size(3) == num_cells);
	this->t4_slope(2, cell) = 0.0;
	this->t4_slope(3, cell) = 0.0;
    }

    return;
}
	
} // end of rtt_imc

//---------------------------------------------------------------------------//
//                              end of Mesh_Operations.cc
//---------------------------------------------------------------------------//
