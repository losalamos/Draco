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

//===========================================================================//
// Constructors and destructors
//===========================================================================//
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

    // Construct an CAR_CU_Mesh int CCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ccsf_i_(long & mesh_index, long & self, long & data, 
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::CCSF<int> > ccsf_i;

	if (idata_size == 0)
	{
	    // Construct a new CAR_CU_Mesh int CCSF class object.
	    ccsf_i = new CAR_CU_Mesh::CCSF<int>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->num_cells(), 
	        "Invalid data passed to construct_mesh_ccsf_i_data_!");

	    vector<int> data_set(mesh->num_cells());

	    for (int cell = 0; cell < mesh->num_cells(); cell++)
	    {
	        data_set[cell] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh int CCSF class object.
	    ccsf_i = new CAR_CU_Mesh::CCSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh CCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::CCSF<int> >::insert(ccsf_i);

    }

    // Destroy a CAR_CU_Mesh int CCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ccsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::CCSF<int> > ccsf_i = 
	    opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	ccsf_i = SP<CAR_CU_Mesh::CCSF<int> >();
	Ensure (!ccsf_i);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::CCSF<int> >::erase(self);
    }

    // Construct an CAR_CU_Mesh double CCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ccsf_d_(long & mesh_index, long & self, double & data,
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;
        SP<CAR_CU_Mesh::CCSF<double> > ccsf_d;

	if (idata_size == 0 )
	{
	    // Construct a new CAR_CU_Mesh double CCSF class object.
	    ccsf_d = new CAR_CU_Mesh::CCSF<double>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->num_cells(), 
	        "Invalid data passed to construct_mesh_ccsf_d_data_!");

	    vector<double> data_set(mesh->num_cells());

	    for (int cell = 0; cell < mesh->num_cells(); cell++)
	    {
	        data_set[cell] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh double CCSF class object.
	    ccsf_d = new CAR_CU_Mesh::CCSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh CCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::CCSF<double> >::insert(ccsf_d);

    }

    // Destroy a CAR_CU_Mesh double CCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ccsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::CCSF<double> > ccsf_d = 
	    opaque_pointers<CAR_CU_Mesh::CCSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	ccsf_d = SP<CAR_CU_Mesh::CCSF<double> >();
	Ensure (!ccsf_d);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::CCSF<double> >::erase(self);
    }

    // Construct a CAR_CU_Mesh int FCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcsf_i_(long & mesh_index, long & self, long & data, 
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::FCSF<int> > fcsf_i;
	
	if (idata_size == 0 )
	{
	    // Construct a new CAR_CU_Mesh int FCSF class object.
	    fcsf_i = new CAR_CU_Mesh::FCSF<int>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->num_face_nodes(), 
	        "Invalid data passed to construct_mesh_fcsf_i_data_!");

	    vector<int> data_set(mesh->num_face_nodes());

	    for (int face = 0; face < mesh->num_face_nodes(); face++)
	    {
	        data_set[face] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh int FCSF class object.
	    fcsf_i = new CAR_CU_Mesh::FCSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::FCSF<int> >::insert(fcsf_i);

    }

    // Destroy a CAR_CU_Mesh int FCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCSF<int> > fcsf_i = 
	    opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	fcsf_i = SP<CAR_CU_Mesh::FCSF<int> >();
	Ensure (!fcsf_i);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCSF<int> >::erase(self);
    }

    // Construct a CAR_CU_Mesh double FCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcsf_d_(long & mesh_index, long & self, double & data,
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;
	SP<CAR_CU_Mesh::FCSF<double> > fcsf_d;

	if (idata_size == 0 )
	{
	    // Construct a new CAR_CU_Mesh double FCSF class object.
	    fcsf_d = new CAR_CU_Mesh::FCSF<double>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->num_face_nodes(), 
	        "Invalid data passed to construct_mesh_fcsf_d_data_!");

	    vector<double> data_set(mesh->num_face_nodes());

	    for (int face = 0; face < mesh->num_face_nodes(); face++)
	    {
	        data_set[face] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh double FCSF class object.
	    fcsf_d = new CAR_CU_Mesh::FCSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::FCSF<double> >::insert(fcsf_d);

    }

    // Destroy a CAR_CU_Mesh double FCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCSF<double> > fcsf_d = 
	    opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	fcsf_d = SP<CAR_CU_Mesh::FCSF<double> >();
	Ensure (!fcsf_d);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCSF<double> >::erase(self);
    }

    // Construct a CAR_CU_Mesh int FCDSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcdsf_i_(long & mesh_index, long & self, long & data,
				 long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::FCDSF<int> >  fcdsf_i;

	if (idata_size == 0)
	{
	    // Construct a new CAR_CU_Mesh int FCDSF class object.
	    fcdsf_i = new CAR_CU_Mesh::FCDSF<int>(mesh);
	}
	else
	{
	    Insist(idata_size ==  mesh->num_cells() * 2 * mesh->get_ndim(), 
	        "Invalid data passed to construct_mesh_fcdsf_i_data_!");

	    vector<vector<int> > data_set(mesh->num_cells());

	    for (int cell = 0; cell < mesh->num_cells(); cell++)
	    {
	        data_set[cell].resize(2 * mesh->get_ndim());
		for (int face = 0; face < 2 * mesh->get_ndim(); face++)
		{
		    data_set[cell][face] = * data_array;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh int FCDSF class object.
	    fcdsf_i = new CAR_CU_Mesh::FCDSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCDSF class object 
	// (self).
	self = opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::insert(fcdsf_i);
    }

    // Destroy a CAR_CU_Mesh int FCDSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcdsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCDSF<int> > fcdsf_i = 
	    opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	fcdsf_i = SP<CAR_CU_Mesh::FCDSF<int> >();
	Ensure (!fcdsf_i);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::erase(self);
    }

    // Construct a CAR_CU_Mesh double FCDSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcdsf_d_(long & mesh_index, long & self, double & data,
				 long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;
	SP<CAR_CU_Mesh::FCDSF<double> > fcdsf_d;

	if (idata_size == 0)
	{
	    // Construct a new CAR_CU_Mesh double FCDSF class object.
	    fcdsf_d = new CAR_CU_Mesh::FCDSF<double>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->num_cells() * 2 * mesh->get_ndim(),
	        "Invalid data passed to construct_mesh_fcdsf_d_data_!");

	    vector<vector<double> > data_set(mesh->num_cells());

	    for (int cell = 0; cell < mesh->num_cells(); cell++)
	    {
	        data_set[cell].resize(2 * mesh->get_ndim());
		for (int face = 0; face < 2 * mesh->get_ndim(); face++)
		{
		    data_set[cell][face] = * data_array;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh double FCDSF class object.
	    fcdsf_d = new CAR_CU_Mesh::FCDSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCDSF class object 
	// (self).
	self = opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::insert(fcdsf_d);
    }

    // Destroy a CAR_CU_Mesh double FCDSF class object from a Fortran 90 
    // program call.
    void destruct_mesh_fcdsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCDSF<double> > fcdsf_d = 
	    opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	fcdsf_d = SP<CAR_CU_Mesh::FCDSF<double> >();
	Ensure (!fcdsf_d);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::erase(self);
    }

//===========================================================================//
// General mesh scalar accessor functions
//===========================================================================//
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

//===========================================================================//
// Layout accessor functions
//===========================================================================//
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

    // Return the face-centered cell node specified by the face.
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

    // Return an array of the nodes that make up a cell, including both the
    // corner nodes and the face-centered nodes.
    void get_mesh_cell_nodes_(long & self, long & cell, long & nodes, 
                              long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_nodes_!");

	vector<int> node_set = mesh->cell_nodes(icell);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of cell nodes (corner + face-centered!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    // Return an array of the face-centered nodes for a cell.
    void get_mesh_cell_face_cen_nodes_(long & self, long & cell, long & nodes,
				       long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell passed in get_mesh_cell_face_centered_nodes_!");

	vector<int> node_set = mesh->cell_face_centered_nodes(icell);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of face-centered cell nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    // Return an array of the corner nodes for a cell.
    void get_mesh_cell_corner_nodes_(long & self, long & cell, long & nodes, 
				     long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_corner_nodes_!");

	vector<int> node_set = mesh->cell_corner_nodes(icell);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of corner cell nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

    //  Return an array of the nodes that comprise a cell face.
    void get_mesh_cell_face_nodes_(long & self, long & cell, long & face,
				   long & nodes, long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_nodes_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_face_nodes_!");

	vector<int> node_set = mesh->cell_face_nodes(icell, iface);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of cell face nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

//===========================================================================//
// Vertex accessor functions
//===========================================================================//
    // Return the entire node vertex array (including both the corner and 
    // face-centered nodes).
    void get_mesh_vertices_(long & self, double & vertices, long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	vector<vector<double> > vertex_set = mesh->get_vertex();

	Insist(ivertex_size == vertex_set.size() * vertex_set[0].size(),
	       "Vertex size error in get_mesh_vertices_!");

	for (int node = 0; node < mesh->num_nodes(); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

    // Return an array containing the vertices for all of the cell corner 
    // nodes.
    void get_mesh_corner_node_vertices_(long & self, double & vertices, 
					long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	vector<vector<double> > vertex_set = mesh->get_vertex();

	Insist(ivertex_size == vertex_set.size() * (vertex_set[0].size() -
	       mesh->num_face_nodes()), 
	       "Vertex size error in get_mesh_corner_node_vertices_!");

	for (int node = 0; node < mesh->num_corner_nodes(); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

    // Return an array containing the vertices for all of the face-centered 
    // nodes.
    void get_mesh_face_cen_node_vertices_(long & self, double & vertices, 
					  long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	vector<vector<double> > vertex_set = mesh->get_vertex();

	Insist(ivertex_size == vertex_set.size() * (vertex_set[0].size() -
	       mesh->num_corner_nodes()), 
	       "Vertex size error in get_mesh_face_centered_node_vertices_!");

	for (int node = mesh->num_corner_nodes(); node < mesh->num_nodes(); 
	     node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

    // Return an array with all of a cell's vertices.
    void get_mesh_cell_vertices_(long & self, long & cell, double & vertices, 
				 long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_vertices_!");

	vector<vector<double> > vertex_set = mesh->get_vertices(icell);

	Insist(ivertex_size == vertex_set.size() * vertex_set[0].size(), 
	       "Vertex size error in get_mesh_cell_vertices_!");

	for (int node = 0; node < pow(2.0, mesh->get_ndim()); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

    // Return an array with all of a cell face's vertices.
    void get_mesh_cell_face_vertices_(long & self, long & cell, long & face, 
				      double & vertices, long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_vertices_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_cell_face_vertices_!");


	vector<vector<double> > vertex_set = mesh->get_vertices(icell, iface);

	Insist(ivertex_size == vertex_set.size() * vertex_set[0].size(), 
	       "Vertex size error in get_mesh_cell_vertices_!");

	for (int node = 0; node < 2.0 * (mesh->get_ndim() - 1); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

    // Return a single node's vertices
    void get_mesh_node_vertices_(long & self, long & node, double & vertices,
				 long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int inode = static_cast<int>(node);
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	Insist(inode > 0 && inode <= mesh->num_nodes(), 
	       "Invalid node number passed to get_mesh_node_vertices_!");
	Insist(ivertex_size == mesh->get_ndim(), 
	       "Vertex size error in get_mesh_node_vertices_!");

	vector<double> vertex_set = mesh->get_vertex(inode);

	for (int dir = 0; dir < mesh->get_ndim(); dir++)
	{
	        * vertex_array = vertex_set[dir];
		++vertex_array;
	}
    }

//===========================================================================//
// Mesh geometry scalar accessor functions
//===========================================================================//
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
	       "Invalid cell number passed to get_mesh_cell_max_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_max_coord_!");

	maximum_value = mesh->max(idir, icell);
    }

    // Return the width in a given direction for the cell in the mesh (self).
    void get_mesh_cell_width_(long & self, long & cell, long & direction, 
				  double & width)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_cellwidth_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_width_!");

	width = mesh->dim(idir, icell);
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

//===========================================================================//
// Mesh field accessor functions
//===========================================================================//
//===========================================================================//
// integer CCSF class objects
//===========================================================================//
    // Return an entire C++ CAR_CU_Mesh integer CCSF class object (self).
    void get_mesh_ccsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<int> ccsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->num_cells(), 
	       "Invalid data passed to get_mesh_ccsf_i_!");

	for (int cell = 0; cell < mesh->num_cells(); cell++)
	{
	    * data_array = static_cast<long>(ccsf_i(cell));
	    ++data_array;
	}

    }

    // Return a cell value from a C++ CAR_CU_Mesh integer CCSF class object
    // (self).
    void get_mesh_ccsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<int> ccsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_ccsf_i_cell_!");

	data = static_cast<long>(ccsf_i(icell));

    }

    // Set an entire C++ CAR_CU_Mesh integer CCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ccsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<int> ccsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->num_cells(), 
	       "Invalid data passed to set_mesh_ccsf_i_!");

	for (int cell = 0; cell < mesh->num_cells(); cell++)
	{
	    ccsf_i(cell) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

    // Set a cell value for a C++ CAR_CU_Mesh integer CCSF class object
    // (self).
    void set_mesh_ccsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<int> ccsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to set_mesh_ccsf_i_cell_!");

	ccsf_i(icell) = static_cast<int>(data);

    }

//===========================================================================//
// double CCSF class objects
//===========================================================================//
    // Return an entire C++ CAR_CU_Mesh integer CCSF class object (self).
    void get_mesh_ccsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<int> ccsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->num_cells(), 
	       "Invalid data passed to get_mesh_ccsf_d_!");

	for (int cell = 0; cell < mesh->num_cells(); cell++)
	{
	    * data_array = ccsf_d(cell);
	    ++data_array;
	}

    }

    // Return a cell value from a C++ CAR_CU_Mesh integer CCSF class object
    // (self).
    void get_mesh_ccsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<int> ccsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to get_mesh_ccsf_d_cell_!");

	data = ccsf_d(icell);

    }

    // Set an entire C++ CAR_CU_Mesh integer CCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ccsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<int> ccsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->num_cells(), 
	       "Invalid data passed to set_mesh_ccsf_d_!");

	for (int cell = 0; cell < mesh->num_cells(); cell++)
	{
	    ccsf_d(cell) = * data_array;
	    ++data_array;
	}
    }

    // Set a cell value for a C++ CAR_CU_Mesh integer CCSF class object
    // (self).
    void set_mesh_ccsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<int> ccsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to set_mesh_ccsf_d_cell__!");

	ccsf_d(icell) = data;

    }

//===========================================================================//
// integer FCSF class objects
//===========================================================================//
    // Return an entire C++ CAR_CU_Mesh integer FCSF class object (self).
    void get_mesh_fcsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<int> fcsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->num_face_nodes(), 
	       "Invalid data passed to get_mesh_fcsf_i_!");

	for (int face = 1; face <= mesh->num_face_nodes(); face++)
	{
	    * data_array = static_cast<long>(fcsf_i(face));
	    ++data_array;
	}

    }

    // Return all of the face values for a cell in a C++ CAR_CU_Mesh integer
    // FCSF class object (self).
    void get_mesh_fcsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> fcsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(),
	       "Invalid data passed to get_mesh_fcsf_i_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = static_cast<long>(fcsf_i(icell, iface));
	    ++data_array;
	}
    }

    // Return a cell face value from a C++ CAR_CU_Mesh integer FCSF class 
    // object (self).
    void get_mesh_fcsf_i_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind, 
				    long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> fcsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  get_mesh_fcsf_i_cell_face_!");

	data = static_cast<long>(fcsf_i(icell, iface));

    }

    // Set an entire C++ CAR_CU_Mesh integer FCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<int> fcsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->num_face_nodes(), 
	       "Invalid data passed to set_mesh_fcsf_i_!");

	for (int face = 1; face <= mesh->num_face_nodes(); face++)
	{
	    fcsf_i(face) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

    // Set all of the face values for a cell in a C++ CAR_CU_Mesh integer FCSF
    // class object (self).
    void set_mesh_fcsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> fcsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data passed to set_mesh_fcsf_i_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    fcsf_i(icell, iface) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

    // Set a cell face value for a C++ CAR_CU_Mesh integer FCSF class object
    // (self).
    void set_mesh_fcsf_i_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
			            long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> fcsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  set_mesh_fcsf_i_cell_face_!");

	fcsf_i(icell, iface) = static_cast<int>(data);

    }

//===========================================================================//
// double FCSF class objects
//===========================================================================//
    // Return an entire C++ CAR_CU_Mesh integer FCSF class object (self).
    void get_mesh_fcsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<double> fcsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->num_face_nodes(), 
	       "Invalid data passed to get_mesh_fcsf_d_!");

	for (int face = 1; face <= mesh->num_face_nodes(); face++)
	{
	    * data_array = fcsf_d(face);
	    ++data_array;
	}

    }

    // Return all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCSF class object (self).
    void get_mesh_fcsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> fcsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data passed to get_mesh_fcsf_d_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = fcsf_d(icell, iface);
	    ++data_array;
	}
    }

    // Return a cell face value for a C++ CAR_CU_Mesh integer FCSF class object
    // (self).
    void get_mesh_fcsf_d_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
				    double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> fcsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  get_mesh_fcsf_d_cell_face_!");

	data = fcsf_d(icell, iface);

    }

    // Set an entire C++ CAR_CU_Mesh integer FCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<double> fcsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->num_face_nodes(), 
	       "Invalid data passed to set_mesh_fcsf_d_!");

	for (int face = 1; face <= mesh->num_face_nodes(); face++)
	{
	    fcsf_d(face) = * data_array;
	    ++data_array;
	}
    }

    // Set all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCSF class object (self).
    void set_mesh_fcsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> fcsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data passed to set_mesh_fcsf_d_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    fcsf_d(icell, iface) = * data_array;
	    ++data_array;
	}
    }

    // Set a cell face value for a C++ CAR_CU_Mesh integer FCSF class object
    // (self).
    void set_mesh_fcsf_d_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
				    double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> fcsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  set_mesh_fcsf_d_cell_face_!");

	fcsf_d(icell, iface) = data;

    }

//===========================================================================//
// integer FCDSF class objects
//===========================================================================//
    // Return an entire C++ CAR_CU_Mesh integer FCDSF class object (self).
    void get_mesh_fcdsf_i_(long & mesh_index, long & self, 
			   long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<int> fcdsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data passed to get_mesh_fcdsf_i_!");

	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        * data_array = static_cast<long>(fcdsf_i(cell, face));
		++data_array;
	    }
	}
    }

    // Return all of the face values for a cell in a C++ CAR_CU_Mesh integer
    // FCDSF class object (self).
    void get_mesh_fcdsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			        long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> fcdsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcdsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(),
	       "Invalid data passed to get_mesh_fcdsf_i_cell_!");

	vector<int> data_set = fcdsf_i(icell);

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = static_cast<long>(data_set[iface]);
	    ++data_array;
	}
    }

    // Return a cell face value from a C++ CAR_CU_Mesh integer FCDSF class 
    // object (self).
    void get_mesh_fcdsf_i_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind, 
				     long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> fcdsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcdsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  get_mesh_fcdsf_i_cell_face_!");

	data = static_cast<long>(fcdsf_i(icell, iface));

    }

    // Set an entire C++ CAR_CU_Mesh integer FCDSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcdsf_i_(long & mesh_index, long & self, 
			   long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<int> fcdsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data passed to set_mesh_fcdsf_i_!");

	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        fcdsf_i(cell, face) = static_cast<int>(* data_array);
		++data_array;
	    }
	}
    }

    // Set all of the face values for a cell in a C++ CAR_CU_Mesh integer FCDSF
    // class object (self).
    void set_mesh_fcdsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			        long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> fcdsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcdsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data passed to set_mesh_fcdsf_i_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
            fcdsf_i(icell, iface) = static_cast<int>(* data_array);
     	    ++data_array;
	}
    }

    // Set a cell face value for a C++ CAR_CU_Mesh integer FCDSF class object
    // (self).
    void set_mesh_fcdsf_i_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind,
			             long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> fcdsf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcdsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  set_mesh_fcdsf_i_cell_face_!");

	fcdsf_i(icell, iface) = static_cast<int>(data);

    }

//===========================================================================//
// double FCDSF class objects
//===========================================================================//
    // Return an entire C++ CAR_CU_Mesh integer FCDSF class object (self).
    void get_mesh_fcdsf_d_(long & mesh_index, long & self, 
			   double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<double> fcdsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size ==  mesh->num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data passed to get_mesh_fcdsf_d_!");

	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        * data_array = fcdsf_d(cell, face);
		++data_array;
	    }
	}
    }

    // Return all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCDSF class object (self).
    void get_mesh_fcdsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			        double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> fcdsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcdsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data passed to get_mesh_fcdsf_d_cell_!");

	vector<double> data_set = fcdsf_d(icell);

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = data_set[iface];
	    ++data_array;
	}
    }

    // Return a cell face value for a C++ CAR_CU_Mesh integer FCDSF class
    // object (self).
    void get_mesh_fcdsf_d_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind,
				     double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> fcdsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  get_mesh_fcdsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  get_mesh_fcdsf_d_cell_face_!");

	data = fcdsf_d(icell, iface);

    }

    // Set an entire C++ CAR_CU_Mesh integer FCDSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcdsf_d_(long & mesh_index, long & self, 
			   double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<double> fcdsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to integer
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size ==  mesh->num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data passed to set_mesh_fcdsf_d_!");

	for (int cell = 1; cell <= mesh->num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        fcdsf_d(cell, face) = * data_array;
		++data_array;
	    }
	}
    }

    // Set all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCDSF class object (self).
    void set_mesh_fcdsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			        double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> fcdsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcdsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data passed to set_mesh_fcdsf_d_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    fcdsf_d(icell, iface) = * data_array;
	    ++data_array;
	}
    }

    // Set a cell face value for a C++ CAR_CU_Mesh integer FCDSF class object
    // (self).
    void set_mesh_fcdsf_d_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind,
				     double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> fcdsf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to integer
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->num_cells(), 
	       "Invalid cell number passed to  set_mesh_fcdsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to  set_mesh_fcdsf_d_cell_face_!");

	fcdsf_d(icell, iface) = data;

    }

} // end extern "C"


} // end namespace rtt_mc

#endif                          // __mc_Shadow_CAR_CU_Mesh_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_CAR_CU_Mesh.cc
//---------------------------------------------------------------------------//
