//----------------------------------*-C++-*----------------------------------//
// Shadow_CAR_CU_Mesh.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
//---------------------------------------------------------------------------//
// @> Shadow_CAR_CU_Mesh interface file
//---------------------------------------------------------------------------//

#ifndef __mc_Shadow_CAR_CU_Mesh_cc__
#define __mc_Shadow_CAR_CU_Mesh_cc__

#include "CAR_CU_Mesh.hh"
#include "Shadow_Opaque_Pointers.hh"
#include "ds++/Assert.hh"
#include <iostream>

//===========================================================================//
// Shadow_CAR_CU_Mesh - 
//
// Purpose : Provides shadow interface functions to the Continuous Adaptive 
// Refinement Cartesion Unstructured Mesh Class for use with Fortran 90 codes.
// Note that the class constructor is not shadowed because the mesh is
// constructed by the Shadow_CAR_CU_Builder class object.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_mc 
{
// draco components
using dsxx::SP;

extern "C" 
{
//---------------------------------------------------------------------------//
// CAR_CU_Mesh F90 to C++ flat interface functions
//---------------------------------------------------------------------------//

    // Destroy a CAR_CU_Mesh class object from a Fortran 90 program call.
    void destruct_car_cu_mesh_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	mesh = SP<CAR_CU_Mesh>();
	Ensure (!mesh);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh>::erase(self);
    }

    // Return the dimension of the mesh (self).
    void get_mesh_dimension_(long & self, long & dimension)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	dimension = mesh->get_ndim();
    }

     // Return the number of cells (ncells) in the mesh (self).
    void get_mesh_num_cells_(long & self, long & ncells)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	ncells = mesh->num_cells();
    }

   // Return the total number of nodes (nnodes) in the mesh (self).
    void get_mesh_num_nodes_(long & self, long & nnodes)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	nnodes = mesh->num_nodes();
    }

    // Return the number of cell-corner nodes (ncnodes) in the mesh (self).
    void get_mesh_num_corner_nodes_(long & self, long & ncnodes)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

        ncnodes = mesh->num_corner_nodes();
    }

    // Return the number of face-centered nodes (nfnodes) in the mesh (self).
    void get_mesh_num_face_nodes_(long & self, long & nfnodes)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	nfnodes = mesh->num_face_nodes();
    }

    // Return the number of cells that are adjacent to this cell face in the 
    // mesh (self).
    void get_mesh_num_adj_(long & self, long & cell, long & face, 
			   long & num_adj)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);	

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_num_adj_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_num_adj_!");

	num_adj = mesh->num_adj(icell, iface);
    }

    // Return the cell that is adjacent to this cell face in the mesh (self).
    void get_mesh_next_cell_(long & self, long & cell, long & face, 
			     long & adj_cell)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_next_cell_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_next_cell_!");

	adj_cell = mesh->next_cell(icell, iface);
    }

    // Return the nth (index) cell that is adjacent to this cell face in the
    // mesh (self). This function is called for cells that are adjacent to 
    // multiple cells.
    void get_mesh_next_specific_cell_(long & self, long & cell, long & face, 
				      long & index, long & adj_cell)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	int iindex = static_cast<int>(index);	

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_next_specific_cell_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_next_specific_cell_!");

	adj_cell = mesh->next_cell(icell, iface, iindex);
    }

    // Return the cell node specified by the index.
    void get_mesh_cell_node_(long & self, long & cell, long & node_index, 
                              long & node)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int ind = static_cast<int>(node_index);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_node_!");
	Insist(ind > 0 && ind <= (2 * mesh->get_ndim() + 
				  pow(2.0,mesh->get_ndim())),
	       "Invalid node index passed to get_mesh_cell_node_!");

	node = mesh->cell_node(icell,ind);
    }

    // Return the face_centered cell node specified by the face.
    void get_mesh_cell_face_cen_node_(long & self, long & cell, long & face, 
                              long & node)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_node_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_face_cen_node_!");

	node = mesh->cell_node(icell,iface);
    }

    // Return the set of nodes that make up a cell, including both the corner
    // nodes and the face-centered nodes.
    void get_mesh_cell_nodes_(long & self, long & cell, long & nodes, 
                              long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_nodes_!");

	vector<int> node_set = mesh->cell_nodes(icell);

	Insist(nodes_size == node_set.size(), 
	       "Illegal number of cell nodes (corner + face-centered!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    // Return the set of corner nodes for a cell.
    void get_mesh_cell_corner_nodes_(long & self, long & cell, long & nodes, 
				     long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_corner_nodes_!");

	vector<int> node_set = mesh->cell_corner_nodes(icell);

	Insist(nodes_size == node_set.size(), 
	       "Illegal number of corner cell nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    // Return the set of face-centered nodes for a cell.
    void get_mesh_cell_face_cen_nodes_(long & self, long & cell, long & nodes,
				       long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell passed in get_mesh_cell_face_centered_nodes_!");

	vector<int> node_set = mesh->cell_face_centered_nodes(icell);

	Insist(nodes_size == node_set.size(), 
	       "Illegal number of face-centered cell nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    // Return the set of face-centered nodes for a cell.
    void get_mesh_cell_face_nodes_(long & self, long & cell, long & face,
				    long & nodes, long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_nodes_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_face_nodes_!");

	vector<int> node_set = mesh->cell_face_nodes(icell, iface);

	Insist(nodes_size == node_set.size(), 
	       "Illegal number of cell face nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    // Return the volume of the cell in the mesh (self).
    void get_mesh_cell_volume_(long & self, long & cell, double & vol)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_volume_!");

	vol = mesh->volume(icell);
    }

    // Return the face area of the cell in the mesh (self).
    void get_mesh_cell_face_area_(long & self, long & cell, long & face,
				  double & area)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_area_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_cell_face_area_!");

	area = mesh->face_area(icell, iface);
    }

    // Return the minimum coordinate value in a given direction for the mesh
    // (self).
    void get_mesh_min_coordinates_(long & self, long & direction, 
				   double & minimum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int idirection = static_cast<int>(direction);

	Insist(idirection > 0 && idirection <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_min_coordinates_!");

	minimum_value = mesh->begin(idirection);
    }

    // Return the maximum coordinate value in a given direction for the mesh
    // (self).
    void get_mesh_max_coordinates_(long & self, long & direction, 
				   double & maximum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int idirection = static_cast<int>(direction);

	Insist(idirection > 0 && idirection <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_max_coordinates_!");

	maximum_value = mesh->end(idirection);
    }

    // Return the minimum coordinate value in a given direction for the cell 
    // in the mesh (self).
    void get_mesh_cell_min_coord_(long & self, long & cell, long & direction, 
				  double & minimum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_min_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_min_coord_!");

	minimum_value = mesh->min(idir, icell);
    }

    // Return the maximum coordinate value in a given direction for the cell 
    // in the mesh (self).
    void get_mesh_cell_max_coord_(long & self, long & cell, long & direction, 
				  double & maximum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cellmax_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_max_coord_!");

	maximum_value = mesh->max(idir, icell);
    }

    // Return the midpoint (i.e., center point) coordinate value in a given
    // direction for a cell in the mesh (self).
    void get_mesh_cell_mid_coord_(long & self, long & cell, long & direction, 
				  double & midpoint_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_mid_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_mid_coord_!");

	midpoint_value = mesh->pos(idir, icell);
    }

    // Return the cell generation level in the mesh (self).
    void get_mesh_cell_generation_(long & self, long & cell, long & generation)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_generation_!");

	generation = mesh->get_generation(icell);
    }

} // end extern "C"


} // end namespace rtt_mc

#endif                          // __mc_Shadow_CAR_CU_Mesh_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Mesh.cc
//---------------------------------------------------------------------------//
