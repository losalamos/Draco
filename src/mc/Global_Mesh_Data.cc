//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Global_Mesh_Data.cc
 * \author Thomas M. Evans
 * \date   Thu Dec  4 17:36:07 2003
 * \brief  Global_Mesh_Data member specializations.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "AMR_Layout.hh"
#include "Layout.hh"
#include "OS_Mesh.hh"
#include "RZWedge_Mesh.hh"
#include "Sphyramid_Mesh.hh"
#include "Global_Mesh_Data.hh"

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// OS_MESH SPECIALIZATION.
//---------------------------------------------------------------------------//
/*!
 * \brief Calculates global data on an OS_Mesh.
 * 
 * For OS_Mesh the spatial extents are:
 * - 0: low x
 * - 1: high x
 * - 2: low y
 * - 3: high y
 * - 4: low z (if 3D)
 * - 5: high z (if 3D)
 * .
 */
template<>
void Global_Mesh_Data<OS_Mesh>::calc_global_mesh_data(const OS_Mesh &mesh)
{
    Require (mesh.num_cells() > 0);

    // calculate spatial extents
    spatial_extents.resize(2 * mesh.get_spatial_dimension());

    // fill them up on this processor and do communication if necessary (DD)
    for (int d = 1, j = 0; d <= mesh.get_spatial_dimension(); d++)
    {
	// get dimensions
	double low  = mesh.begin(d);
	double high = mesh.end(d);

	// find the min and max over all processors
	if (topology->get_parallel_scheme() == "DD")
	{
	    rtt_c4::global_min(low);
	    rtt_c4::global_max(high);
	}
	else if (topology->get_parallel_scheme() != "replication")
	{
	    throw rtt_dsxx::assertion(
		"Invalid parallel scheme in Global_Mesh_Data!");
	}

	// assign to the extents data
	spatial_extents[j++] = low;
	spatial_extents[j++] = high;

	Check (low < high);
    }   
}

//---------------------------------------------------------------------------//
// RZWEDGE_MESH SPECIALIZATION.
//---------------------------------------------------------------------------//
/*!
 * \brief Calculates global data on an RZWedge_Mesh.
 * - 0: low x
 * - 1: high x
 * - 2: low z
 * - 3: high z
 * .
 */
template<>
void Global_Mesh_Data<RZWedge_Mesh>::calc_global_mesh_data(
    const RZWedge_Mesh &mesh)
{
    Require (mesh.num_cells() > 0);

    // resize spatial extents
    spatial_extents.resize(4);

    // pick a cell and search orthogonally in the low x, high x, low z, high
    // z directions until a problem or processor boundary is hit
    bool hit_boundary = false;
    int  next_cell    = 0;
    
    // guessed cell; this may not be optimal depending upon the mesh ordering
    int cell = 1;
   
    // get the AMR_Layout
    const AMR_Layout &layout = mesh.get_Layout();

    // the low x direction is zero by definition in the RZWedge_Mesh
    spatial_extents[0] = 0.0;

    // search the high x direction
    while (!hit_boundary)
    {
	// get the next cell in the high x direction
	next_cell = layout(cell, 2, 1);

	// check for boundary
	//   reflection        : next_cell = cell
	//   vacuum            : next_cell = 0
	//   processor boundary: next_cell < 0
	if (next_cell == cell || next_cell <= 0)
	{
	    // assign high x position
	    spatial_extents[1] = mesh.get_high_x(cell);
	    Check (spatial_extents[1] > 0.0);

	    // end search
	    hit_boundary = true;
	}
	else
	{
	    // update to next cell
	    cell = next_cell;
	    Check (next_cell > 0);
	}
    }

    // search the low z direction
    hit_boundary = false;
    cell         = 1;
    while (!hit_boundary)
    {
	// get the next cell in the low z direction
	next_cell = layout(cell, 5, 1);

	// check for boundary
	//   reflection        : next_cell = cell
	//   vacuum            : next_cell = 0
	//   processor boundary: next_cell < 0
	if (next_cell == cell || next_cell <= 0)
	{
	    // assign low z position
	    spatial_extents[2] = mesh.get_low_z(cell);

	    // end search
	    hit_boundary = true;
	}
	else
	{
	    // update to next cell
	    cell = next_cell;
	    Check (next_cell > 0);
	}
    }

    // search the high z direction
    hit_boundary = false;
    cell         = 1;
    while (!hit_boundary)
    {
	// get the next cell in the high z direction
	next_cell = layout(cell, 6, 1);

	// check for boundary
	//   reflection        : next_cell = cell
	//   vacuum            : next_cell = 0
	//   processor boundary: next_cell < 0
	if (next_cell == cell || next_cell <= 0)
	{
	    // assign low z position
	    spatial_extents[3] = mesh.get_high_z(cell);

	    // end search
	    hit_boundary = true;
	}
	else
	{
	    // update to next cell
	    cell = next_cell;
	    Check (next_cell > 0);
	}
    }
    
    // find the min and max over all processors
    if (topology->get_parallel_scheme() == "DD")
    {
	// get global minimum for low z
	rtt_c4::global_min(spatial_extents[2]);

	// get global maximum for high x and high z
	rtt_c4::global_max(spatial_extents[1]);
	rtt_c4::global_max(spatial_extents[3]);
    }

    Check (spatial_extents[1] > 0.0);
    Check (spatial_extents[2] < spatial_extents[3]);
    Check (spatial_extents[0] == 0.0);
}

//---------------------------------------------------------------------------//
// SPHYRAMID_MESH SPECIALIZATION.
//---------------------------------------------------------------------------//
/*!
 * \brief Calculates global data on an Sphyramid_Mesh.
 * - 0: low x
 * - 1: high x
 * .
 */
template<>
void Global_Mesh_Data<Sphyramid_Mesh>::calc_global_mesh_data(
    const Sphyramid_Mesh &mesh)
{
    Require (mesh.num_cells() > 0);

    // resize spatial extents
    spatial_extents.resize(2);

    // pick a cell and search orthogonally in the high x directions until a
    // problem or processor boundary is hit
    bool hit_boundary = false;
    int  next_cell    = 0;
    
    // guessed cell; this may not be optimal depending upon the mesh ordering
    int cell = 1;
   
    // get the AMR_Layout
    const Layout &layout = mesh.get_Layout();

    // the low x direction is zero by definition in the Sphyramid_Mesh
    spatial_extents[0] = 0.0;

    // search the high x direction
    while (!hit_boundary)
    {
	// get the next cell in the high x direction
	next_cell = layout(cell, 2);

	// check for boundary
	//   reflection        : next_cell = cell
	//   vacuum            : next_cell = 0
	//   processor boundary: next_cell < 0
	if (next_cell == cell || next_cell <= 0)
	{
	    // assign high x position
	    spatial_extents[1] = mesh.get_high_x(cell);
	    Check (spatial_extents[1] > 0.0);

	    // end search
	    hit_boundary = true;
	}
	else
	{
	    // update to next cell
	    cell = next_cell;
	    Check (next_cell > 0);
	}
    }
    
    // find the max over all processors
    if (topology->get_parallel_scheme() == "DD")
    {
	// get global maximum for high x
	rtt_c4::global_max(spatial_extents[1]);
    }

    Check (spatial_extents[0] == 0.0);
    Check (spatial_extents[1] > 0.0);
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                 end of Global_Mesh_Data.cc
//---------------------------------------------------------------------------//
