//----------------------------------*-C++-*----------------------------------//
// Shadow_Mesh.cc
// B.T. Adams (bta@lanl.gov)
// 27 Sept 99
/*! 
 * \file   amr_mesh/Shadow_Mesh.cc
 * \author B.T. Adams
 * \date   Mon 27 Sep 10:33:26 1999
 * \brief  Provides the C++ side of the shadow object interface functions to
 *         the Continuous Adaptive Refinement Cartesion Unstructured Mesh 
 *         class for use with Fortran 90 codes. The complimentary Fortran 90
 *         shadow object interface functions that reference the functions 
 *         herein are provided in Shadow_Mesh.f90. Note that the class
 *         constructor is not shadowed because the mesh is constructed 
 *         directly by the Shadow_Builder class object.  An example 
 *         code is also provide via link to illustrate the usage of all of 
 *         the shadow object interface functions to the amr_mesh package from
 *         a Fortran 90 code.
 *
 *\sa Mark G. Gray, Randy M. Roberts, and Tom Evans, Scientific Programming,
 *   "Shadow-Object Interface Between Fortran 95 and C++", March-April 1999.
 */
//---------------------------------------------------------------------------//
// @> Shadow_Mesh interface file
//---------------------------------------------------------------------------//

#ifndef __amr_Shadow_Mesh_cc__
#define __amr_Shadow_Mesh_cc__

#include "Mesh.hh"
#include "ds++/opaquePointers.hh"
#include "ds++/Assert.hh"
#include <iostream>

//===========================================================================//
// Shadow_Mesh - 
//
// Purpose : Provides shadow object interface functions to the Continuous 
//           Adaptive Refinement Cartesion Unstructured Mesh class for use
//           with Fortran 90 codes. Note that the class constructor is not 
//           shadowed because the mesh is constructed directly by the 
//           Shadow_Builder class object.
//
// revision history:
// -----------------
//  0) original
// 
//===========================================================================//

namespace rtt_amr 
{
// draco components
using rtt_dsxx::SP;
using rtt_dsxx::opaque_pointers;

extern "C" 
{
//---------------------------------------------------------------------------//
// CAR_CU_Mesh F90 to C++ flat interface functions
//---------------------------------------------------------------------------//

//===========================================================================//
// Constructors and destructors
//===========================================================================//
/*!
 * \brief Shadow object that destroys the CAR_CU_Mesh class object that
 *        is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 */
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

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int cell-centered scalar field (CCSF) object sized to the number 
 *        of cells in the mesh from a Fortran 90 program call. The object will
 *        be initialized if data_size > 0.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int CCSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             int CCSF object (returned).
 * \param data Initialization int data for the CCSF (optional)
 * \param data_size Size of the CCSF initialization data vector (optional, but
 *                  must equal the number of cells in the mesh if provided).
 */
    // Construct an CAR_CU_Mesh int CCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ccsf_i_(long & mesh_index, long & self, long & data, 
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::CCSF<int> > CCSF_int;

	if (idata_size == 0)
	{
	    // Construct a new CAR_CU_Mesh int CCSF class object.
	    CCSF_int = new CAR_CU_Mesh::CCSF<int>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_cells(), 
	        "Invalid data size passed to construct_mesh_ccsf_i_!");

	    vector<int> data_set(mesh->get_num_cells());

	    for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	    {
	        data_set[cell] = static_cast<int>(* data_array) ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh int CCSF class object.
	    CCSF_int = new CAR_CU_Mesh::CCSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh CCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::CCSF<int> >::insert(CCSF_int);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        int cell-centered scalar field (CCSF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class
 *             int CCSF object.
 */
    // Destroy a CAR_CU_Mesh int CCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ccsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::CCSF<int> > CCSF_int = 
	    opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	CCSF_int = SP<CAR_CU_Mesh::CCSF<int> >();
	Ensure (!CCSF_int);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::CCSF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double cell-centered scalar field (CCSF) object sized to the number
 *        of cells in the mesh from a Fortran 90 program call. The object will
 *        be initialized if data_size > 0.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double CCSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             double CCSF object (returned).
 * \param data Initialization double data for the CCSF (optional)
 * \param data_size Size of the CCSF initialization data vector (optional, but
 *                  must equal the number of cells in the mesh if provided).
 */
    // Construct an CAR_CU_Mesh double CCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ccsf_d_(long & mesh_index, long & self, double & data,
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;
        SP<CAR_CU_Mesh::CCSF<double> > CCSF_double;

	if (idata_size == 0 )
	{
	    // Construct a new CAR_CU_Mesh double CCSF class object.
	    CCSF_double = new CAR_CU_Mesh::CCSF<double>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_cells(), 
	        "Invalid data size passed to construct_mesh_ccsf_d_!");

	    vector<double> data_set(mesh->get_num_cells());

	    for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	    {
	        data_set[cell] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh double CCSF class object.
	    CCSF_double = new CAR_CU_Mesh::CCSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh CCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::CCSF<double> >::insert(CCSF_double);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class 
 *        double cell-centered scalar field (CCSF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class CCSF
 *             object.
 */
    // Destroy a CAR_CU_Mesh double CCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ccsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::CCSF<double> > CCSF_double = 
	    opaque_pointers<CAR_CU_Mesh::CCSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	CCSF_double = SP<CAR_CU_Mesh::CCSF<double> >();
	Ensure (!CCSF_double);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::CCSF<double> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int cell-centered vector field (CCVF) object from a Fortran 90 
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the first array index is specified by lead_indx. The 
 *        dimension of the second array index is fixed to be the same as the
 *        number of cells in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int CCVF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             int CCVF object (returned).
 * \param data Initialization int data for the CCVF (optional)
 * \param data_size Size of the CCVF initialization data vector (optional,
 *                  but must equal the number of cells in the mesh times the 
 *                  leading index if provided).
 * \param lead_indx Size of the CCVF initialization data vector leading index.
 */
    // Construct an CAR_CU_Mesh int CCVF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0. The
    // dimension of the first array index is specified by lead_index and 
    // defaults to the problem geometry dimension. The dimension of the second
    // array index is fixed to be the same as the number of cells. Note that
    // this is the exact opposite of what occurs on the Fortran side, where 
    // the array dimensions are assumed to be ncells x lead_indx (and the 1D 
    // array data is passed to this routine with that assumption).
    void construct_mesh_ccvf_i_(long & mesh_index, long & self, long & data, 
				long & data_size, long & lead_indx)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ilead_indx = static_cast<int>(lead_indx);
	long * data_array = & data;
	SP<CAR_CU_Mesh::CCVF<int> > CCVF_int;

	if (idata_size == 0)
	{
	    if (ilead_indx == mesh->get_ndim())
	        // Construct a new CAR_CU_Mesh int CCVF class object with the
	        // default of the size of the vector leading index equal to 
	        // that of the problem geometry.
	        CCVF_int = new CAR_CU_Mesh::CCVF<int>(mesh);
	    else
	        // Construct a new CAR_CU_Mesh int CCVF class object with the
	        // size of the vector leading index arbitrary
	        CCVF_int = new CAR_CU_Mesh::CCVF<int>(mesh, ilead_indx);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_cells() * ilead_indx, 
	           "Invalid data size passed to construct_mesh_ccvf_i_!");

	    vector<vector<int> > data_set(ilead_indx);

	    for (int dim = 0; dim < ilead_indx; dim++)
	        data_set[dim].resize(mesh->get_num_cells());

	    for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	    {	    
	        for (int dim = 0; dim < ilead_indx; dim++)
		{
		    data_set[dim][cell] = static_cast<int>(* data_array) ;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh int CCVF class object.
	    CCVF_int = new CAR_CU_Mesh::CCVF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh CCVF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::CCVF<int> >::insert(CCVF_int);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class 
 *        int cell-centered vector field (CCVF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class int
 *             CCVF object.
 */
    // Destroy a CAR_CU_Mesh int CCVF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ccvf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::CCVF<int> > CCVF_int = 
	    opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	CCVF_int = SP<CAR_CU_Mesh::CCVF<int> >();
	Ensure (!CCVF_int);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::CCVF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double cell-centered vector field (CCVF) object from a Fortran 90 
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the first array index is specified by lead_indx. The 
 *        dimension of the second array index is fixed to be the same as the
 *        number of cells in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double CCVF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             double CCVF object (returned).
 * \param data Initialization double data for the CCVF (optional)
 * \param data_size Size of the CCVF initialization data vector (optional,
 *                  but must equal the number of cells in the mesh times the 
 *                  leading index if provided).
 * \param lead_indx Size of the CCVF initialization data vector leading index.
 */
    // Construct an CAR_CU_Mesh double CCVF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0. The
    // dimension of the first array index is specified by lead_index and 
    // defaults to the problem geometry dimension. The dimension of the second
    // array index is fixed to be the same as the number of cells. Note that
    // this is the exact opposite of what occurs on the Fortran side, where 
    // the array dimensions are assumed to be ncells x lead_index (and the 1D 
    // array data is passed to this routine with that assumption).
    void construct_mesh_ccvf_d_(long & mesh_index, long & self, double & data,
				long & data_size, long & lead_indx)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ilead_indx = static_cast<int>(lead_indx);
	double * data_array = & data;
        SP<CAR_CU_Mesh::CCVF<double> > CCVF_double;

	if (idata_size == 0 )
	{
	    if (ilead_indx == mesh->get_ndim())
	        // Construct a new CAR_CU_Mesh int CCVF class object with the
	        // default of the size of the vector leading index equal to 
	        // that of the problem geometry.
	        CCVF_double = new CAR_CU_Mesh::CCVF<double>(mesh);
	    else
	        // Construct a new CAR_CU_Mesh int CCVF class object with the
	        // size of the vector leading index arbitrary
	        CCVF_double = new CAR_CU_Mesh::CCVF<double>(mesh, ilead_indx);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_cells() * ilead_indx, 
	        "Invalid data size passed to construct_mesh_ccvf_d_!");

	    vector<vector<double> > data_set(ilead_indx);

	    for (int dim = 0; dim < ilead_indx; dim++)
	        data_set[dim].resize(mesh->get_num_cells());

	    for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	    {
	        for (int dim = 0; dim < ilead_indx; dim++)
		{
		    data_set[dim][cell] = * data_array ;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh double CCVF class object.
	    CCVF_double = new CAR_CU_Mesh::CCVF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh CCVF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::CCVF<double> >::insert(CCVF_double);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        double CCVF object from a Fortran 90 program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 */
    // Destroy a CAR_CU_Mesh double CCVF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ccvf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::CCVF<double> > CCVF_double = 
	    opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	CCVF_double = SP<CAR_CU_Mesh::CCVF<double> >();
	Ensure (!CCVF_double);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::CCVF<double> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int face-centered scalar field (FCSF) object sized to the number of
 *        unique cell faces in the mesh from a Fortran 90 program call. The 
 *        object will be initialized if data_size > 0.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             int FCSF object (returned).
 * \param data Initialization int data for the FCSF (optional)
 * \param data_size Size of the FCSF initialization data vector (optional,
 *                  but must equal the number of cell faces in the mesh if 
 *                  provided).
 */
    // Construct a CAR_CU_Mesh int FCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcsf_i_(long & mesh_index, long & self, long & data, 
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::FCSF<int> > FCSF_int;
	
	if (idata_size == 0 )
	{
	    // Construct a new CAR_CU_Mesh int FCSF class object.
	    FCSF_int = new CAR_CU_Mesh::FCSF<int>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_face_nodes(), 
	        "Invalid data size passed to construct_mesh_fcsf_i_!");

	    vector<int> data_set(mesh->get_num_face_nodes());

	    for (int face = 0; face < mesh->get_num_face_nodes(); face++)
	    {
	        data_set[face] = static_cast<int>(* data_array) ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh int FCSF class object.
	    FCSF_int = new CAR_CU_Mesh::FCSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::FCSF<int> >::insert(FCSF_int);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        int  face-centered scalar field (FCSF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class int
 *             FCSF object.
 */
    // Destroy a CAR_CU_Mesh int FCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCSF<int> > FCSF_int = 
	    opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	FCSF_int = SP<CAR_CU_Mesh::FCSF<int> >();
	Ensure (!FCSF_int);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCSF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double face-centered scalar field (FCSF) object sized to the number
 *        of unique cell faces in the mesh from a Fortran 90 program call. The 
 *        object will be initialized if data_size > 0.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class
 *             double FCSF object (returned).
 * \param data Initialization double data for the FCSF (optional)
 * \param data_size Size of the FCSF initialization data vector (optional,
 *                  but must equal the number of cell faces in the mesh if 
 *                  provided).
 */
    // Construct a CAR_CU_Mesh double FCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcsf_d_(long & mesh_index, long & self, double & data,
				long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;
	SP<CAR_CU_Mesh::FCSF<double> > FCSF_double;

	if (idata_size == 0 )
	{
	    // Construct a new CAR_CU_Mesh double FCSF class object.
	    FCSF_double = new CAR_CU_Mesh::FCSF<double>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_face_nodes(), 
	        "Invalid data size passed to construct_mesh_fcsf_d_!");

	    vector<double> data_set(mesh->get_num_face_nodes());

	    for (int face = 0; face < mesh->get_num_face_nodes(); face++)
	    {
	        data_set[face] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh double FCSF class object.
	    FCSF_double = new CAR_CU_Mesh::FCSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::FCSF<double> >::insert(FCSF_double);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        double face-centered scalar field (FCSF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 */
    // Destroy a CAR_CU_Mesh double FCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCSF<double> > FCSF_double = 
	    opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	FCSF_double = SP<CAR_CU_Mesh::FCSF<double> >();
	Ensure (!FCSF_double);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCSF<double> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int face-centered discontinuous scalar field (FCDSF) class object
 *        sized to the number of discontinuous cell faces in the mesh from
 *        a Fortran 90 program call. The object will be initialized if 
 *        data_size > 0.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCDSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object (returned).
 * \param data Initialization int data for the FCDSF (optional)
 * \param data_size Size of the FCDSF initialization data vector (optional,
 *                  but must equal the number of discontinuous cell faces in 
 *                  the mesh if provided).
 */
    // Construct a CAR_CU_Mesh int FCDSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcdsf_i_(long & mesh_index, long & self, long & data,
				 long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::FCDSF<int> >  FCDSF_int;

	if (idata_size == 0)
	{
	    // Construct a new CAR_CU_Mesh int FCDSF class object.
	    FCDSF_int = new CAR_CU_Mesh::FCDSF<int>(mesh);
	}
	else
	{
	    Insist(idata_size ==  mesh->get_num_cells() * 2 * mesh->get_ndim(), 
	        "Invalid data size passed to construct_mesh_fcdsf_i_!");

	    vector<vector<int> > data_set(mesh->get_num_cells());

	    for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	    {
	        data_set[cell].resize(2 * mesh->get_ndim());
		for (int face = 0; face < 2 * mesh->get_ndim(); face++)
		{
		    data_set[cell][face] = static_cast<int>(* data_array);
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh int FCDSF class object.
	    FCDSF_int = new CAR_CU_Mesh::FCDSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCDSF class object 
	// (self).
	self = opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::insert(FCDSF_int);
    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        int face-centered discontinuous scalar field (FCDSF) object from a
 *        Fortran 90 program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class int
 *             FCDSF object.
 */
    // Destroy a CAR_CU_Mesh int FCDSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcdsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCDSF<int> > FCDSF_int = 
	    opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	FCDSF_int = SP<CAR_CU_Mesh::FCDSF<int> >();
	Ensure (!FCDSF_int);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double face-centered discontinuous scalar field (FCDSF) class object
 *        sized to the number of discontinuous cell faces in the mesh from
 *        a Fortran 90 program call. The object will be initialized if 
 *        data_size > 0.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCDSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class
 *             double FCDSF object (returned).
 * \param data Initialization double data for the FCDSF (optional)
 * \param data_size Size of the FCDSF initialization data vector (optional,
 *                  but must equal the number of discontinuous cell faces in 
 *                  the mesh if provided).
 */    
    // Construct a CAR_CU_Mesh double FCDSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcdsf_d_(long & mesh_index, long & self, double & data,
				 long & data_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;
	SP<CAR_CU_Mesh::FCDSF<double> > FCDSF_double;

	if (idata_size == 0)
	{
	    // Construct a new CAR_CU_Mesh double FCDSF class object.
	    FCDSF_double = new CAR_CU_Mesh::FCDSF<double>(mesh);
	}
	else
	{
	    Insist(idata_size == mesh->get_num_cells() * 2 * mesh->get_ndim(),
	        "Invalid data size passed to construct_mesh_fcdsf_d_!");

	    vector<vector<double> > data_set(mesh->get_num_cells());

	    for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	    {
	        data_set[cell].resize(2 * mesh->get_ndim());
		for (int face = 0; face < 2 * mesh->get_ndim(); face++)
		{
		    data_set[cell][face] = * data_array;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh double FCDSF class object.
	    FCDSF_double = new CAR_CU_Mesh::FCDSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCDSF class object 
	// (self).
	self = opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::insert(FCDSF_double);
    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        double face-centered discontinuous scalar field (FCDSF) object from
 *        a Fortran 90 program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 */
    // Destroy a CAR_CU_Mesh double FCDSF class object from a Fortran 90 
    // program call.
    void destruct_mesh_fcdsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCDSF<double> > FCDSF_double = 
	    opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	FCDSF_double = SP<CAR_CU_Mesh::FCDSF<double> >();
	Ensure (!FCDSF_double);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int face-centered vector field (FCVF) object from a Fortran 90 
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the second array index is specified by vec_size. The
 *        dimension of the first array index is fixed to be the same as the
 *        number of unique cell faces in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCVF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class
 *             int FCVF object (returned).
 * \param data Initialization int data for the FCVF (optional)
 * \param data_size Size of the FCVF initialization data vector (optional,
 *                  but must equal the number of cell faces in the mesh times
 *                  the trailing index if provided).
 * \param vec_size Size of the FCVF initialization data vector trailing index.
 */
    // Construct a CAR_CU_Mesh int FCVF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcvf_i_(long & mesh_index, long & self, long & data, 
				long & data_size, long & vec_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ivec_size = static_cast<int>(vec_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::FCVF<int> > fcvf_i;

	if (idata_size == 0 )
	{
	    if (ivec_size == mesh->get_ndim())
	    {
	        // Construct a new CAR_CU_Mesh int FCVF class object with the
	        // dimension of the second array element being the same as 
	        // that of the problem geometry.
	        fcvf_i = new CAR_CU_Mesh::FCVF<int>(mesh);
	    }
	    else
	    {
	        // Construct a new CAR_CU_Mesh int FCVF class object with the
	        // dimension of the second array element specified by vec_size.
	        fcvf_i = new CAR_CU_Mesh::FCVF<int>(mesh, vec_size);
	    }
	}
	else
	{
	    Insist(idata_size == mesh->get_num_face_nodes() * ivec_size, 
	           "Invalid data size passed to construct_mesh_fcvf_i_!");

	    vector<vector<int> > data_set(mesh->get_num_face_nodes());

	    for (int face = 0; face < mesh->get_num_face_nodes(); face++)
	    {
	        data_set[face].resize(vec_size);
	        for (int dim = 0; dim < vec_size; dim++)
		{
		    data_set[face][dim] = static_cast<int>(* data_array) ;
		    ++data_array;
		}
	    }
	    // Construct a new CAR_CU_Mesh int FCVF class object.
	    fcvf_i = new CAR_CU_Mesh::FCVF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCVF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::FCVF<int> >::insert(fcvf_i);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        int  face-centered vector field (FCVF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class int
 *             FCVF object.
 */
    // Destroy a CAR_CU_Mesh int FCVF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcvf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCVF<int> > fcvf_i = 
	    opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	fcvf_i = SP<CAR_CU_Mesh::FCVF<int> >();
	Ensure (!fcvf_i);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCVF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double face-centered vector field (FCVF) object from a Fortran 90 
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the second array index is specified by vec_size. The
 *        dimension of the first array index is fixed to be the same as the
 *        number of unique cell faces in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCVF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class
 *             double FCVF object (returned).
 * \param data Initialization double data for the FCVF (optional)
 * \param data_size Size of the FCVF initialization data vector (optional,
 *                  but must equal the number of cell faces in the mesh times
 *                  the trailing index if provided).
 * \param vec_size Size of the FCVF initialization data vector trailing index.
 */
    // Construct a CAR_CU_Mesh double FCVF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_fcvf_d_(long & mesh_index, long & self, double & data,
				long & data_size, long & vec_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ivec_size = static_cast<int>(vec_size);
	double * data_array = & data;
	SP<CAR_CU_Mesh::FCVF<double> > fcvf_d;

	if (idata_size == 0 )
	{
	    if (ivec_size == mesh->get_ndim())
	    {
	        // Construct a new CAR_CU_Mesh double FCVF class object with 
	        // the dimension of the second array element being the same as 
	        // that of the problem geometry.
	        fcvf_d = new CAR_CU_Mesh::FCVF<double>(mesh);
	    }
	    else
	    {
	        // Construct a new CAR_CU_Mesh double FCVF class object with 
	        // the dimension of the second array element  specified by 
	        // vec_size
	        fcvf_d = new CAR_CU_Mesh::FCVF<double>(mesh, ivec_size);
	    }
	}
	else
	{
	    Insist(idata_size == mesh->get_num_face_nodes() * ivec_size, 
	           "Invalid data size passed to construct_mesh_fcvf_d_!");

	   vector< vector<double> > data_set(mesh->get_num_face_nodes());

	    for (int face = 0; face < mesh->get_num_face_nodes(); face++)
	    {
	        data_set[face].resize(vec_size);
	        for (int dim = 0; dim < vec_size; dim++)
		{
		    data_set[face][dim] = * data_array ;
		    ++data_array;
		}
	    }
	    // Construct a new CAR_CU_Mesh double FCVF class object.
	    fcvf_d = new CAR_CU_Mesh::FCVF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh FCVF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::FCVF<double> >::insert(fcvf_d);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        double face-centered vector field (FCVF) object from a Fortran 90
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 */
    // Destroy a CAR_CU_Mesh double FCVF class object from a Fortran 90 program
    // call.
    void destruct_mesh_fcvf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::FCVF<double> > fcvf_d = 
	    opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	fcvf_d = SP<CAR_CU_Mesh::FCVF<double> >();
	Ensure (!fcvf_d);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::FCVF<double> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int node-centered scalar field (NCSF) object from a Fortran 90
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the array is specified by vec_size and must equal 
 *        either the total number of cell-corner and face-centered nodes in
 *        the mesh or the number of cell-corner nodes in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int NCSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             int NCSF object (returned).
 * \param data Initialization int data for the NCSF (optional)
 * \param data_size Size of the NCSF initialization data vector (optional, but
 *                  must equal either the total number of cell-corner and 
 *                  face-centered nodes in the mesh or the number of 
 *                  cell- corner nodes in the mesh if provided).
 * \param vec_size  Size of the NCSF data vector (must equal either the total
 *                  number of cell-corner and face-centered nodes in the mesh
 *                  or the number of cell-corner nodes in the mesh).
 */
    // Construct an CAR_CU_Mesh int NCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ncsf_i_(long & mesh_index, long & self, long & data, 
				long & data_size, long & vec_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ivec_size = static_cast<int>(vec_size);
	long * data_array = & data;
	SP<CAR_CU_Mesh::NCSF<int> > NCSF_int;

	if (idata_size == 0 && ivec_size == mesh->get_num_nodes())
	{
	    // Construct a new unitialized CAR_CU_Mesh int NCSF class object 
	    // for all of the nodes
	    NCSF_int = new CAR_CU_Mesh::NCSF<int>(mesh);
	}
	else if (idata_size == 0 && ivec_size != mesh->get_num_nodes())
	{
	    Insist(ivec_size == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncsf_i_!");

	    // Construct a new unitialized CAR_CU_Mesh int NCSF class object 
	    // for the cell-corner nodes.
	    NCSF_int = new CAR_CU_Mesh::NCSF<int>(mesh, ivec_size);
	}
	else
	{
	    Insist(ivec_size == mesh->get_num_nodes() ||
		   ivec_size == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncsf_i_!");

	    Insist(idata_size == ivec_size, 
	        "Invalid data size passed to construct_mesh_ncsf_i_!");

	    vector<int> data_set(ivec_size);

	    for (int node = 0; node < ivec_size; node++)
	    {
	        data_set[node] = static_cast<int>(* data_array) ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh int NCSF class object.
	    NCSF_int = new CAR_CU_Mesh::NCSF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh NCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::NCSF<int> >::insert(NCSF_int);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        int node-centered scalar field (NCSF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class int 
 *             NCSF object.
 */
    // Destroy a CAR_CU_Mesh int NCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ncsf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::NCSF<int> > NCSF_int = 
	    opaque_pointers<CAR_CU_Mesh::NCSF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	NCSF_int = SP<CAR_CU_Mesh::NCSF<int> >();
	Ensure (!NCSF_int);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::NCSF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double node-centered scalar field (NCSF) object from a Fortran 90
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the array is specified by vec_size and must equal 
 *        either the total number of cell-corner and face-centered nodes in
 *        the mesh or the number of cell-corner nodes in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double NCSF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             double NCSF object (returned).
 * \param data Initialization double data for the NCSF (optional)
 * \param data_size Size of the NCSF initialization data vector (optional, but
 *                  must equal either the total number of cell-corner and 
 *                  face-centered nodes in the mesh or the number of 
 *                  cell-corner nodes in the mesh if provided).
 * \param vec_size  Size of the NCSF data vector (must equal either the total
 *                  number of cell-corner and face-centered nodes in the mesh
 *                  or the number of cell-corner nodes in the mesh).
 */
    // Construct an CAR_CU_Mesh double NCSF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ncsf_d_(long & mesh_index, long & self, double & data,
				long & data_size, long & vec_size)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ivec_size = static_cast<int>(vec_size);
	double * data_array = & data;
        SP<CAR_CU_Mesh::NCSF<double> > NCSF_double;

	if (idata_size == 0 && ivec_size == mesh->get_num_nodes())
	{
	    // Construct a new unitialized CAR_CU_Mesh double NCSF class 
	    // object for all of the nodes
	    NCSF_double = new CAR_CU_Mesh::NCSF<double>(mesh);
	}
	else if (idata_size == 0 && ivec_size != mesh->get_num_nodes())
	{
	    Insist(ivec_size == mesh->get_num_corner_nodes(), 
	           "Invalid data size passed to construct_mesh_ncsf_d_!");

	    // Construct a new unitialized CAR_CU_Mesh double NCSF class 
	    // object for the cell-corner nodes
	    NCSF_double = new CAR_CU_Mesh::NCSF<double>(mesh, ivec_size);
	}
	else
	{
	    Insist(ivec_size == mesh->get_num_nodes() ||
		   ivec_size == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncsf_d_!");

	    Insist(idata_size == ivec_size, 
	        "Invalid data size passed to construct_mesh_ncsf_d_!");

	    vector<double> data_set(ivec_size);

	    for (int node = 0; node < ivec_size; node++)
	    {
	        data_set[node] = * data_array ;
		++data_array;
	    }

	    // Construct a new CAR_CU_Mesh double NCSF class object.
	    NCSF_double = new CAR_CU_Mesh::NCSF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh NCSF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::NCSF<double> >::insert(NCSF_double);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        double node-centered scalar field (NCSF) class object from a Fortran
 *        90 program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCSF object.
 */
    // Destroy a CAR_CU_Mesh double NCSF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ncsf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::NCSF<double> > NCSF_double = 
	    opaque_pointers<CAR_CU_Mesh::NCSF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	NCSF_double = SP<CAR_CU_Mesh::NCSF<double> >();
	Ensure (!NCSF_double);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::NCSF<double> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        int node-centered vector field (NCVF) object from a Fortran 90 
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the second array index is specified by vec_size_2. The
 *        dimension of the first array index is specified by vec_size_1 and
 *        must equal either the total number of cell-corner and face-centered
 *        nodes in the mesh or the number of cell-corner nodes in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int NCVF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             int NCVF object (returned).
 * \param data Initialization int data for the NCVF (optional)
 * \param data_size Size of the NCVF initialization data vector (optional,
 *                  but must equal the number of cell faces in the mesh times
 *                  the trailing index if provided).
 * \param vec_size_1 Size of the NCVF initialization data vector leading index
 *                   (must equal either the total number of cell-corner and 
 *                   face-centered nodes in the mesh or the number of 
 *                   cell-corner nodes in the mesh).
 * \param vec_size_2 Size of the NCVF initialization data vector trailing 
 *                   index.
 */
    // Construct an CAR_CU_Mesh int NCVF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ncvf_i_(long & mesh_index, long & self, long & data, 
				long & data_size, long & vec_size_1,
				long & vec_size_2)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ivec_size_1 = static_cast<int>(vec_size_1);
	int ivec_size_2 = static_cast<int>(vec_size_2);
	long * data_array = & data;
	SP<CAR_CU_Mesh::NCVF<int> > NCVF_int;

	if (idata_size == 0)
	{
	    if (ivec_size_1 == mesh->get_num_nodes() && 
		ivec_size_2 == mesh->get_ndim())
	    {
	        // Construct a new unitialized CAR_CU_Mesh int NCVF class 
	        // object for all of the nodes with the size of the second 
	        // array index defaulted to be equal to the problem geometry 
	        // dimension.
	        NCVF_int = new CAR_CU_Mesh::NCVF<int>(mesh);
	    }
	    else if (ivec_size_1 != mesh->get_num_nodes() && 
		     ivec_size_2 == mesh->get_ndim())
	    {
	        Insist(ivec_size_1 == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncvf_i_!");

	        // Construct a new unitialized CAR_CU_Mesh int NCVF class 
	        // object for the cell-corner nodes with the size of the 
		// second array index defaulted to be equal to the problem 
		// geometry dimension.
	        NCVF_int = new CAR_CU_Mesh::NCVF<int>(mesh, ivec_size_1);
	    }
	    else
	    {
	        Insist(ivec_size_1 == mesh->get_num_nodes() ||
		       ivec_size_1 == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncvf_i_!");

	        // Construct a new unitialized CAR_CU_Mesh int NCVF class 
	        // object for the cell-corner nodes with the size of the 
		// second array index arbitrary.
	        NCVF_int = new CAR_CU_Mesh::NCVF<int>(mesh, ivec_size_1, 
						    ivec_size_2);
	    }
	}
	else
	{
	    Insist(ivec_size_1 == mesh->get_num_nodes() ||
		   ivec_size_1 == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncvf_i_!");

	    Insist(idata_size == ivec_size_1 * ivec_size_2, 
	        "Invalid data size passed to construct_mesh_ncvf_i_!");

	    vector<vector<int> > data_set(ivec_size_1);

	    for (int node = 0; node < ivec_size_1; node++)
	    {
	        data_set[node].resize(ivec_size_2);
	        for (int dim = 0; dim < ivec_size_2; dim++)
		{
		    data_set[node][dim] = static_cast<int>(* data_array) ;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh int NCVF class object.
	    NCVF_int = new CAR_CU_Mesh::NCVF<int>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh NCVF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::NCVF<int> >::insert(NCVF_int);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        int node-centered vector field (NCVF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class int
 *             NCVF object.
 */
    // Destroy a CAR_CU_Mesh int NCVF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ncvf_i_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::NCVF<int> > NCVF_int = 
	    opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	NCVF_int = SP<CAR_CU_Mesh::NCVF<int> >();
	Ensure (!NCVF_int);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::NCVF<int> >::erase(self);
    }

/*!
 * \brief Shadow object that constructs a CAR_CU_Mesh nested mesh field class
 *        double node-centered vector field (NCVF) object from a Fortran 90 
 *        program call. The object will be initialized if data_size > 0. The
 *        dimension of the second array index is specified by vec_size_2. The
 *        dimension of the first array index is specified by vec_size_1 and
 *        must equal either the total number of cell-corner and face-centered
 *        nodes in the mesh or the number of cell-corner nodes in the mesh.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double NCVF class object.
 * \param self Opaque pointer to the new CAR_CU_Mesh nested mesh field class 
 *             double NCVF object (returned).
 * \param data Initialization double data for the NCVF (optional)
 * \param data_size Size of the NCVF initialization data vector (optional,
 *                  but must equal the number of cell faces in the mesh times
 *                  the trailing index if provided).
 * \param vec_size_1 Size of the NCVF initialization data vector leading index
 *                   (must equal either the total number of cell-corner and 
 *                   face-centered nodes in the mesh or the number of 
 *                   cell-corner nodes in the mesh).
 * \param vec_size_2 Size of the NCVF initialization data vector trailing 
 *                   index.
 */
    // Construct an CAR_CU_Mesh double NCVF class object from a Fortran 90 
    // program call. The object will be initialized if data_size > 0.
    void construct_mesh_ncvf_d_(long & mesh_index, long & self, double & data,
				long & data_size, long & vec_size_1,
				long & vec_size_2)
    {
	// Get the address of the CAR_CU_Mesh class object (mesh).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	int ivec_size_1 = static_cast<int>(vec_size_1);
	int ivec_size_2 = static_cast<int>(vec_size_2);
	double * data_array = & data;
        SP<CAR_CU_Mesh::NCVF<double> > NCVF_double;

	if (idata_size == 0)
	{
	    if (ivec_size_1 == mesh->get_num_nodes() && 
		ivec_size_2 == mesh->get_ndim())
	    {
	        // Construct a new unitialized CAR_CU_Mesh int NCVF class 
	        // object for all of the nodes with the size of the second 
	        // array index defaulted to be equal to the problem geometry 
	        // dimension.
	        NCVF_double = new CAR_CU_Mesh::NCVF<double>(mesh);
	    }
	    else if (ivec_size_1 != mesh->get_num_nodes() && 
		     ivec_size_2 == mesh->get_ndim())
	    {
	        Insist(ivec_size_1 == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncvf_d_!");

	        // Construct a new unitialized CAR_CU_Mesh int NCVF class 
	        // object for the cell-corner nodes with the size of the 
		// second array index defaulted to be equal to the problem
		// geometry dimension.
	        NCVF_double = new CAR_CU_Mesh::NCVF<double>(mesh, ivec_size_1);
	    }
	    else
	    {
	        Insist(ivec_size_1 == mesh->get_num_nodes() ||
		       ivec_size_1 == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncvf_d_!");

	        // Construct a new unitialized CAR_CU_Mesh int NCVF class 
	        // object for the cell-corner nodes with the size of the
		// second array index arbitrary.
	        NCVF_double = new CAR_CU_Mesh::NCVF<double>(mesh, ivec_size_1, 
						       ivec_size_2);
	    }
	}
	else
	{
	    Insist(ivec_size_1 == mesh->get_num_nodes() ||
		   ivec_size_1 == mesh->get_num_corner_nodes(), 
	           "Invalid vector size passed to construct_mesh_ncvf_d_!");

	    Insist(idata_size == ivec_size_1 * ivec_size_2, 
	        "Invalid data size passed to construct_mesh_ncvf_d_!");

	    vector<vector<double> > data_set(ivec_size_1);

	    for (int node = 0; node < ivec_size_1; node++)
	    {
	        data_set[node].resize(ivec_size_2);
	        for (int dim = 0; dim < ivec_size_2; dim++)
		{
		    data_set[node][dim] = static_cast<double>(* data_array) ;
		    ++data_array;
		}
	    }

	    // Construct a new CAR_CU_Mesh int NCVF class object.
	    NCVF_double = new CAR_CU_Mesh::NCVF<double>(mesh, data_set);
	}

	// return the map key for the new CAR_CU_Mesh NCVF class object (self).
	self = opaque_pointers<CAR_CU_Mesh::NCVF<double> >::insert(NCVF_double);

    }

/*!
 * \brief Shadow object that destroys a CAR_CU_Mesh nested mesh field class
 *        double node-centered vector field (NCVF) object from a Fortran 90 
 *        program call.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class
 *             double NCVF object.
 */
    // Destroy a CAR_CU_Mesh double NCVF class object from a Fortran 90 program
    // call.
    void destruct_mesh_ncvf_d_(long & self)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh::NCVF<double> > NCVF_double = 
	    opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);

	// destroy the CAR_CU_Mesh class object by assigning this SP to 
        // a null SP
	NCVF_double = SP<CAR_CU_Mesh::NCVF<double> >();
	Ensure (!NCVF_double);

	// remove the opaque pointer to the CAR_CU_Mesh class object.
	opaque_pointers<CAR_CU_Mesh::NCVF<double> >::erase(self);
    }

//===========================================================================//
// General mesh scalar accessor functions
//===========================================================================//
/*!
 * \brief Shadow object that returns the number of spatial dimensions for 
 *        the CAR_CU_Mesh class object that is referenced by the specified 
 *        opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param dimension The number of spatial dimensions (returned).
 */
    // Return the dimension of the mesh (self).
    void get_mesh_dimension_(long & self, long & dimension)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	dimension = mesh->get_ndim();
    }

/*!
 * \brief Shadow object that returns the number of cells in the CAR_CU_Mesh
 *        class object that is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param ncells The number of cells (returned).
 */
     // Return the number of cells (ncells) in the mesh (self).
    void get_mesh_num_cells_(long & self, long & ncells)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	ncells = mesh->get_num_cells();
    }

/*!
 * \brief Shadow object that returns the total number of nodes (i.e., both
 *        cell-corner and face-centered nodes) in the CAR_CU_Mesh class 
 *        object that is referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param nnodes The number of nodes (returned).
 */
   // Return the total number of nodes (nnodes) in the mesh (self).
    void get_mesh_num_nodes_(long & self, long & nnodes)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	nnodes = mesh->get_num_nodes();
    }

/*!
 * \brief Shadow object that returns the number of cell-corner nodes in the
 *        CAR_CU_Mesh class object that is referenced by the specified opaque
 *        pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param ncnodes The number of cell-corner nodes (returned).
 */
    // Return the number of cell-corner nodes (ncnodes) in the mesh (self).
    void get_mesh_num_corner_nodes_(long & self, long & ncnodes)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

        ncnodes = mesh->get_num_corner_nodes();
    }

/*!
 * \brief Shadow object that returns the number of face-centered nodes in the
 *        CAR_CU_Mesh class object that is referenced by the specified opaque
 *        pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param nfnodes The number of face-centered nodes (returned).
 */
    // Return the number of face-centered nodes (nfnodes) in the mesh (self).
    void get_mesh_num_face_nodes_(long & self, long & nfnodes)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);

	nfnodes = mesh->get_num_face_nodes();
    }

//===========================================================================//
// Layout accessor functions
//===========================================================================//
/*!
 * \brief Shadow object that returns the number of cells that are adjacent 
 *        to the specified cell face in the CAR_CU_Mesh class object that is
 *        referenced by the specified opaque pointer.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param face The cell face number.
 * \param num_adj The number of cells that are adjacent to the specified cell
 *        face (returned).
 */
    // Return the number of cells that are adjacent to this cell face in the 
    // mesh (self).
    void get_mesh_num_adj_(long & self, long & cell, long & face, 
			   long & num_adj)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);	

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_num_adj_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_num_adj_!");

	num_adj = mesh->get_num_adj_cells(icell, iface);
    }

/*!
 * \brief Shadow object that returns the number of the cell that is adjacent 
 *        to the specified cell face in the CAR_CU_Mesh class object that is
 *        referenced by the specified opaque pointer. This function is called
 *        for cells that are adjacent to a single cell.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param face The cell face number.
 * \param adj_cell The adjacent cell number (returned).
 */
    // Return the cell that is adjacent to this cell face in the mesh (self).
    void get_mesh_next_cell_(long & self, long & cell, long & face, 
			     long & adj_cell)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_next_cell_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_next_cell_!");

	adj_cell = mesh->get_next_cell(icell, iface);
    }

/*!
 * \brief Shadow object that returns the number of the cell with the specified
 *        index that is adjacent to the specified cell face in the CAR_CU_Mesh
 *        class object that is referenced by the specified opaque pointer. 
 *        This function is called for cells that are adjacent to multiple 
 *        cells.
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param face The cell face number.
 * \param index The adjacent cell index number.
 * \param adj_cell The adjacent cell number (returned).
 */
    // Return the nth (index) cell that is adjacent to this cell face in the
    // mesh (self). This function is called for cells that are adjacent to 
    // multiple cells.
    void get_mesh_next_specific_cell_(long & self, long & cell, long & face, 
				      long & index, long & adj_cell)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	int iindex = static_cast<int>(index);	

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_next_specific_cell_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_next_specific_cell_!");

	adj_cell = mesh->get_next_cell(icell, iface, iindex);
    }

/*!
 * \brief Shadow object that returns the node number for the specified cell 
 *        node index in the CAR_CU_Mesh class object that is referenced by
 *        the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param node_index The cell node index number.
 * \param node The node number (returned).
 */
    // Return the cell node specified by the index.
    void get_mesh_cell_node_(long & self, long & cell, long & node_index, 
                             long & node)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int ind = static_cast<int>(node_index);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_node_!");
	Insist(ind > 0 && ind <= (2 * mesh->get_ndim() + 
				  pow(2.0,mesh->get_ndim())),
	       "Invalid node index passed to get_mesh_cell_node_!");

	node = mesh->get_cell_node(icell,ind);
    }

/*!
 * \brief Shadow object that returns the face-centered node number for the 
 *        specified cell face in the CAR_CU_Mesh class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param face The cell face number.
 * \param node The face-centered node number (returned).
 */
    // Return the face-centered cell node specified by the face.
    void get_mesh_cell_face_cen_node_(long & self, long & cell, long & face, 
				      long & node)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_node_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_face_cen_node_!");

	node = mesh->get_cell_node(icell,iface);
    }

/*!
 * \brief Shadow object that returns all of the nodes (i.e., both the 
 *        cell-corner and the face-centered nodes) that make up the specified
 *        cell in the CAR_CU_Mesh class object that is referenced by the
 *        specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param nodes The cell nodes (returned).
 * \param nodes_size The number of cell nodes (must be equal to the number
 *                   of cell-corner and face-centered nodes for the cell type).
 */
    // Return an array of the nodes that make up a cell, including both the
    // cell-corner nodes and the face-centered nodes.
    void get_mesh_cell_nodes_(long & self, long & cell, long & nodes, 
                              long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_nodes_!");

	vector<int> node_set = mesh->get_cell_nodes(icell);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of cell nodes (corner + face-centered!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

/*!
 * \brief Shadow object that returns all of the face-centered nodes that make
 *        up the specified cell in the CAR_CU_Mesh class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param nodes The cell face-centered nodes (returned).
 * \param nodes_size The number of cell face-centered nodes (must be equal 
 *                   to the number of face-centered nodes for the cell type).
 */
    // Return an array of the face-centered nodes for a cell.
    void get_mesh_cell_face_cen_nodes_(long & self, long & cell, long & nodes,
				       long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell passed in get_mesh_cell_face_centered_nodes_!");

	vector<int> node_set = mesh->get_cell_face_centered_nodes(icell);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of face-centered cell nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

/*!
 * \brief Shadow object that returns all of the cell-corner nodes that make
 *        up the specified cell in the CAR_CU_Mesh class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param nodes The cell-corner nodes (returned).
 * \param nodes_size The number of cell-corner nodes (must be equal 
 *                   to the number of cell-corner nodes for the cell type).
 */
    // Return an array of the cell-corner nodes for a cell.
    void get_mesh_cell_corner_nodes_(long & self, long & cell, long & nodes, 
				     long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_corner_nodes_!");

	vector<int> node_set = mesh->get_cell_corner_nodes(icell);

	Insist(inodes_size == node_set.size(), 
	       "Invalid number of corner cell nodes!");

	for (int node_index = 0; node_index < node_set.size(); node_index++)
	{
	    * nodes_array = node_set[node_index];
	    ++nodes_array;
	}
    }

/*!
 * \brief Shadow object that returns all of the face-corner nodes that make
 *        up the specified cell face in the CAR_CU_Mesh class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param face The cell face number.
 * \param nodes The cell face-corner nodes (returned).
 * \param nodes_size The number of cell face-corner nodes (must be equal to
 *                   the number of face-corner nodes for the cell face type).
 */
    //  Return an array of the nodes that comprise a cell face.
    void get_mesh_cell_face_nodes_(long & self, long & cell, long & face,
				   long & nodes, long & nodes_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	int inodes_size = static_cast<int>(nodes_size);
	long * nodes_array = & nodes;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_nodes_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_face_nodes_!");

	vector<int> node_set = mesh->get_cell_face_nodes(icell, iface);

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
/*!
 * \brief Shadow object that returns all of the node coordinates (i.e., each 
 *        spatial direction for both the cell-corner and the face-centered 
 *        nodes) in the CAR_CU_Mesh class object that is referenced by the
 *        specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param vertices The node coordinates (returned).
 * \param vertex_size The size of the node vertex return array (must be equal
 *                    to the the number of spatial dimensions times the total 
 *                    number of nodes).
 */
    // Return the entire node vertex array (including both the cell-corner and 
    // face-centered nodes).
    void get_mesh_vertices_(long & self, double & vertices, long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	vector<vector<double> > vertex_set = mesh->get_node_coords();

	Insist(ivertex_size == vertex_set.size() * vertex_set[0].size(),
	       "Vertex size error in get_mesh_vertices_!");

	for (int node = 0; node < mesh->get_num_nodes(); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the cell-corner node coordinates
 *        (i.e., each spatial direction for the cell-corner nodes) in the 
 *        CAR_CU_Mesh class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param vertices The cell-corner node coordinates (returned).
 * \param vertex_size The size of the cell-corner node vertex return array 
 *                    (must be equal to the the number of spatial dimensions
 *                    times the number of cell-corner nodes).
 */
    // Return an array containing the vertices for all of the cell-corner 
    // nodes.
    void get_mesh_corner_node_vertices_(long & self, double & vertices, 
					long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	vector<vector<double> > vertex_set = mesh->get_node_coords();

	Insist(ivertex_size == vertex_set.size() * (vertex_set[0].size() -
	       mesh->get_num_face_nodes()), 
	       "Vertex size error in get_mesh_corner_node_vertices_!");

	for (int node = 0; node < mesh->get_num_corner_nodes(); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the face-centered node coordinates
 *        (i.e., each spatial direction for the face-centered nodes) in the 
 *        CAR_CU_Mesh class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param vertices The face-centered node coordinates (returned).
 * \param vertex_size The size of the face-centered node vertex return array 
 *                    (must be equal to the the number of spatial dimensions
 *                    times the number of face-centered nodes).
 */
    // Return an array containing the vertices for all of the face-centered 
    // nodes.
    void get_mesh_face_cen_node_vertices_(long & self, double & vertices, 
					  long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	vector<vector<double> > vertex_set = mesh->get_node_coords();

	Insist(ivertex_size == vertex_set.size() * (vertex_set[0].size() -
	       mesh->get_num_corner_nodes()), 
	       "Vertex size error in get_mesh_face_centered_node_vertices_!");

	for (int node = mesh->get_num_corner_nodes(); 
	     node < mesh->get_num_nodes(); node++)
	{
	    for (int dir = 0; dir < mesh->get_ndim(); dir++)
	    {
	        * vertex_array = vertex_set[dir][node];
		++vertex_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the node coordinates (i.e., each 
 *        spatial direction for both the cell-corner and the face-centered 
 *        nodes) for the specified cell in the CAR_CU_Mesh class object that
 *        is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param vertices The cell node coordinates (returned).
 * \param vertex_size The size of the cell node vertex return array (must be
 *                    equal to the the number of spatial dimensions times the
 *                    total number of cell nodes).
 */
    // Return an array with all of a cell's vertices.
    void get_mesh_cell_vertices_(long & self, long & cell, double & vertices, 
				 long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_vertices_!");

	vector<vector<double> > vertex_set = 
	    mesh->get_cell_nodes_coords(icell);

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

/*!
 * \brief Shadow object that returns all of the node coordinates (i.e., each 
 *        spatial direction for the cell face-corner nodes) for the specified
 *        cell face in the CAR_CU_Mesh class object that is referenced by the
 *        specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param cell The cell face number.
 * \param vertices The cell node coordinates (returned).
 * \param vertex_size The size of the cell node vertex return array (must be
 *                    equal to the the number of spatial dimensions times the
 *                    total number of cell face-corner nodes).
 */
    // Return an array with all of a cell face's vertices.
    void get_mesh_cell_face_vertices_(long & self, long & cell, long & face, 
				      double & vertices, long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_vertices_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_cell_face_vertices_!");


	vector<vector<double> > vertex_set = 
	    mesh->get_cell_nodes_coords(icell, iface);

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

/*!
 * \brief Shadow object that returns all of the node coordinates (i.e., each 
 *        spatial direction) for the specified node in the CAR_CU_Mesh class 
 *        object that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param node The node number.
 * \param vertices The node coordinates (returned).
 * \param vertex_size The size of the node vertex return array (must be equal
 *                    to the the number of spatial dimensions).
 */
    // Return a single node's vertices
    void get_mesh_node_vertices_(long & self, long & node, double & vertices,
				 long & vertex_size)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node);
	int ivertex_size = static_cast<int>(vertex_size);
	double * vertex_array = & vertices;

	Insist(inode > 0 && inode <= mesh->get_num_nodes(), 
	       "Invalid node number passed to get_mesh_node_vertices_!");
	Insist(ivertex_size == mesh->get_ndim(), 
	       "Vertex size error in get_mesh_node_vertices_!");

	vector<double> vertex_set = mesh->get_node_coords(inode);

	for (int dir = 0; dir < mesh->get_ndim(); dir++)
	{
	        * vertex_array = vertex_set[dir];
		++vertex_array;
	}
    }

//===========================================================================//
// Mesh geometry scalar accessor functions
//===========================================================================//
/*!
 * \brief Shadow object that returns the volume of the specified cell in 
 *        the CAR_CU_Mesh class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param vol The cell volume (returned).
 */
    // Return the volume of the cell in the mesh (self).
    void get_mesh_cell_volume_(long & self, long & cell, double & vol)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_volume_!");

	vol = mesh->get_cell_volume(icell);
    }

/*!
 * \brief Shadow object that returns the area of the specified cell face in 
 *        the CAR_CU_Mesh class object that is referenced by the specified 
 *        opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param face The cell face number.
 * \param area The cell face area (returned).
 */
    // Return the face area of the cell in the mesh (self).
    void get_mesh_cell_face_area_(long & self, long & cell, long & face,
				  double & area)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int iface = static_cast<int>(face);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_face_area_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_cell_face_area_!");

	area = mesh->get_cell_face_area(icell, iface);
    }

/*!
 * \brief Shadow object that returns the minimum coordinate value along the
 *        specified direction for the CAR_CU_Mesh class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param direction Coordinate direction (x=1, y=2, z =3).
 * \param minimum_value Minimum coordinate value (returned).
 */
    // Return the minimum coordinate value in a given direction for the mesh
    // (self).
    void get_mesh_min_coordinates_(long & self, long & direction, 
				   double & minimum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int idirection = static_cast<int>(direction);

	Insist(idirection > 0 && idirection <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_min_coordinates_!");

	minimum_value = mesh->get_mesh_min_coord(idirection);
    }

/*!
 * \brief Shadow object that returns the maximum coordinate value along the
 *        specified direction for the CAR_CU_Mesh class object that is 
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param direction Coordinate direction (x=1, y=2, z =3).
 * \param maximum_value Maximum coordinate value (returned).
 */
    // Return the maximum coordinate value in a given direction for the mesh
    // (self).
    void get_mesh_max_coordinates_(long & self, long & direction, 
				   double & maximum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int idirection = static_cast<int>(direction);

	Insist(idirection > 0 && idirection <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_max_coordinates_!");

	maximum_value = mesh->get_mesh_max_coord(idirection);
    }

/*!
 * \brief Shadow object that returns the minimum coordinate value along the
 *        specified direction within the specified cell in the CAR_CU_Mesh
 *        class object that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param direction Coordinate direction (x=1, y=2, z =3).
 * \param minimum_value Minimum coordinate value (returned).
 */
    // Return the minimum coordinate value in a given direction for the cell 
    // in the mesh (self).
    void get_mesh_cell_min_coord_(long & self, long & cell, long & direction, 
				  double & minimum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_min_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_min_coord_!");

	minimum_value = mesh->get_cell_min_coord(idir, icell);
    }

/*!
 * \brief Shadow object that returns the mid-point coordinate value along the
 *        specified direction within the specified cell in the CAR_CU_Mesh
 *        class object that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param direction Coordinate direction (x=1, y=2, z =3).
 * \param midpoint_value Mid-point coordinate value (returned).
 */
    // Return the midpoint (i.e., center point) coordinate value in a given
    // direction for a cell in the mesh (self).
    void get_mesh_cell_mid_coord_(long & self, long & cell, long & direction, 
				  double & midpoint_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_mid_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_mid_coord_!");

	midpoint_value = mesh->get_cell_center_coord(idir, icell);
    }

/*!
 * \brief Shadow object that returns the maximum coordinate value along the
 *        specified direction within the specified cell in the CAR_CU_Mesh
 *        class object that is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param direction Coordinate direction (x=1, y=2, z =3).
 * \param maximum_value Maximum coordinate value (returned).
 */
    // Return the maximum coordinate value in a given direction for the cell 
    // in the mesh (self).
    void get_mesh_cell_max_coord_(long & self, long & cell, long & direction, 
				  double & maximum_value)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_max_coord_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_max_coord_!");

	maximum_value = mesh->get_cell_max_coord(idir, icell);
    }

/*!
 * \brief Shadow object that returns the width of the specified cell along
 *        the specified direction within the CAR_CU_Mesh class object that
 *        is referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param direction Coordinate direction (x=1, y=2, z =3).
 * \param width The cell width in the specified direction (returned).
 */
    // Return the width in a given direction for the cell in the mesh (self).
    void get_mesh_cell_width_(long & self, long & cell, long & direction, 
				  double & width)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);
	int idir = static_cast<int>(direction);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cellwidth_!");
	Insist(idir > 0 && idir <= mesh->get_ndim(), 
	       "Invalid direction passed to get_mesh_cell_width_!");

	width = mesh->get_cell_width(idir, icell);
    }

/*!
 * \brief Shadow object that returns the generation (i.e., refinement) level
 *        of the specified cell within the CAR_CU_Mesh class object that is
 *        referenced by the specified opaque pointer. 
 * \param self Opaque pointer to the CAR_CU_Mesh class object.
 * \param cell The cell number.
 * \param generation The cell generation level (returned).
 */
    // Return the cell generation level in the mesh (self).
    void get_mesh_cell_generation_(long & self, long & cell, long & generation)
    {
	// Get the address of the CAR_CU_Mesh class object (self).
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_cell_generation_!");

	generation = mesh->get_generation(icell);
    }

//===========================================================================//
// Mesh field accessor functions
//===========================================================================//
//===========================================================================//
// int CCSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int cell-centered scalar field (CCSF) object that
 *        is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCSF object.
 * \param data CCSF int data values (returned).
 * \param data_size Size of the CCSF returned data vector (must equal the 
 *                  number of cells in the mesh).
 */
    // Return an entire C++ CAR_CU_Mesh int CCSF class object (self).
    void get_mesh_ccsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<int> CCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_cells(), 
	       "Invalid data size passed to get_mesh_ccsf_i_!");

	for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	{
	    * data_array = static_cast<long>(CCSF_int(cell));
	    ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns the specified cell value from the 
 *        CAR_CU_Mesh nested mesh field class int cell-centered scalar field
 *       (CCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCSF object.
 * \param cell_ind The cell number. 
 * \param data CCSF int cell data value (returned).
 */
    // Return a cell value from a C++ CAR_CU_Mesh int CCSF class object
    // (self).
    void get_mesh_ccsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<int> CCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_ccsf_i_cell_!");

	data = static_cast<long>(CCSF_int(icell));

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int cell-centered scalar field (CCSF) object 
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCSF object.
 * \param data CCSF int data values (supplied).
 * \param data_size Size of the CCSF supplied data vector (must equal the 
 *                  number of cells in the mesh).
 */
    // Set an entire C++ CAR_CU_Mesh int CCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ccsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<int> CCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_cells(), 
	       "Invalid data size passed to set_mesh_ccsf_i_!");

	for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	{
	    CCSF_int(cell) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell value for the 
 *        CAR_CU_Mesh nested mesh field class int cell-centered scalar field
 *        (CCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCSF object.
 * \param cell_ind The cell number. 
 * \param data CCSF cell int data value (supplied).
 */
    // Set a cell value for a C++ CAR_CU_Mesh int CCSF class object
    // (self).
    void set_mesh_ccsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<int> CCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_ccsf_i_cell_!");

	CCSF_int(icell) = static_cast<int>(data);

    }

//===========================================================================//
// double CCSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double cell-centered scalar field (CCSF) object 
 *        that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCSF object.
 * \param data CCSF double data values (returned).
 * \param data_size Size of the CCSF returned data vector (must equal the 
 *                  number of cells in the mesh).
 */
    // Return an entire C++ CAR_CU_Mesh double CCSF class object (self).
    void get_mesh_ccsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<double> CCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->get_num_cells(), 
	       "Invalid data size passed to get_mesh_ccsf_d_!");

	for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	{
	    * data_array = CCSF_double(cell);
	    ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns the specified cell value from the 
 *        CAR_CU_Mesh nested mesh field class double cell-centered scalar 
 *        field (CCSF) object that is referenced by the specified opaque
 *        pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCSF object.
 * \param cell_ind The cell number. 
 * \param data CCSF double cell data value (returned).
 */
    // Return a cell value from a C++ CAR_CU_Mesh double CCSF class object
    // (self).
    void get_mesh_ccsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<double> CCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_ccsf_d_cell_!");

	data = CCSF_double(icell);

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double cell-centered scalar field (CCSF) object 
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCSF object.
 * \param data CCSF double data values (supplied).
 * \param data_size Size of the CCSF supplied data vector (must equal the 
 *                  number of cells in the mesh).
 */
    // Set an entire C++ CAR_CU_Mesh double CCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ccsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCSF<double> CCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->get_num_cells(), 
	       "Invalid data size passed to set_mesh_ccsf_d_!");

	for (int cell = 0; cell < mesh->get_num_cells(); cell++)
	{
	    CCSF_double(cell) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell value for the 
 *        CAR_CU_Mesh nested mesh field class double cell-centered scalar 
 *        field (CCSF) object that is referenced by the specified opaque
 *        pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double CCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCSF object.
 * \param cell_ind The cell number. 
 * \param data CCSF double cell data value (supplied).
 */
    // Set a cell value for a C++ CAR_CU_Mesh double CCSF class object
    // (self).
    void set_mesh_ccsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCSF<double> CCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_ccsf_d_cell__!");

	CCSF_double(icell) = data;

    }

//===========================================================================//
// int CCVF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int cell-centered vector field (CCVF) object that
 *        is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCVF object.
 * \param data CCVF int data values (returned)
 * \param data_size Size of the CCVF returned data vector (must equal the 
 *                  number of cells in the mesh times the leading index size).
 */
    // Return an entire C++ CAR_CU_Mesh int CCVF class object (self) - 
    // works for both the arbitrary leading index size and the default with 
    // the size of the leading index equal to that of the problem geometry 
    // dimension.
    void get_mesh_ccvf_i_(long & mesh_index, long & self, long & data, 
			  long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCVF<int> CCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_cells() * CCVF_int.get_size(), 
	       "Invalid data size passed to get_mesh_ccvf_i_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int dim = 1; dim <= CCVF_int.get_size(); dim++)
	    {
	        * data_array = static_cast<long>(CCVF_int(dim, cell));
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the leading index values for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class int
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCVF object.
 * \param cell_ind The cell number.
 * \param data CCVF int cell data values (returned).
 * \param data_size Size of the CCVF returned data vector (must equal the 
 *                  size of the CCVF leading index).
 */
    // Return all of the dim values for a cell in a C++ CAR_CU_Mesh integer
    // CCVF class object (self) - works for both the arbitrary leading index
    // size and the default with the size of the leading index equal to that
    // of the problem geometry dimension.
    void get_mesh_ccvf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<int> CCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_ccvf_i_cell_!");
	Insist(idata_size == CCVF_int.get_size(), 
	       "Invalid data size passed to get_mesh_ccvf_i_cell_!");

	for (int idim = 1; idim <= idata_size; idim++)
	{
	    * data_array = static_cast<long>(CCVF_int(idim, icell));
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified leading index value for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class int 
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCVF object.
 * \param cell_ind The cell number.
 * \param dim_ind The leading index number.
 * \param data CCVF int leading index, cell data value (returned).
 */
    // Return a cell dim value from a C++ CAR_CU_Mesh int CCVF class 
    // object (self) - works for both the arbitrary leading index size 
    // and the default with the size of the leading index equal to that
    // of the problem geometry dimension.
    void get_mesh_ccvf_i_cell_dim_(long & mesh_ind, long & self, 
				   long & cell_ind, long & dim_ind, 
				   long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<int> CCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_ccvf_i_cell_dim_!");
	Insist(idim > 0 && idim <= CCVF_int.get_size(), 
	       "Invalid dim number passed to get_mesh_ccvf_i_cell_dim_!");

	data = static_cast<long>(CCVF_int(idim, icell));

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int cell-centered vector field (CCVF) object that
 *        is referenced by the specified opaque pointers. This can also be
 *        done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCVF object.
 * \param data CCVF int data values (supplied).
 * \param data_size Size of the CCVF supplied data vector (must equal the 
 *                  number of cells in the mesh times the leading index size).
 */
    // Set an entire C++ CAR_CU_Mesh int CCVF class object (self) (can 
    // also be done at initialization using the constructor) - works for 
    // both the arbitrary leading index size and the default with the size 
    // of the leading index equal to that of the problem geometry dimension.
    void set_mesh_ccvf_i_(long & mesh_index, long & self, long & data, 
			  long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCVF<int> CCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_cells() * CCVF_int.get_size(), 
	       "Invalid data size passed to set_mesh_ccvf_i_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int dim = 1; dim <= CCVF_int.get_size(); dim++)
	    {
	        CCVF_int(dim, cell) = static_cast<int>(* data_array);
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that sets all of the leading index values for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class int 
 *        cell-centered vector field (CCVF) object that is referenced by 
 *        the specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCVF object.
 * \param cell_ind The cell number.
 * \param data CCVF int cell data values (supplied).
 * \param data_size Size of the CCVF supplied data vector (must equal the 
 *                  size of the CCVF leading index).
 */
    // Set all of the dim values for a cell in a C++ CAR_CU_Mesh int CCVF
    // class object (self) - works for both the arbitrary leading index size 
    // and the default with the size of the leading index equal to that of the
    // problem geometry dimension.
    void set_mesh_ccvf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<int> CCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_ccvf_i_cell_!");
	Insist(idata_size == CCVF_int.get_size(), 
	       "Invalid data size passed to set_mesh_ccvf_i_cell_!");

	for (int idim = 1; idim <= CCVF_int.get_size(); idim++)
	{
            CCVF_int(idim, icell) = static_cast<int>(* data_array);
     	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified leading index value for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class int 
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int CCVF object.
 * \param cell_ind The cell number.
 * \param dim_ind The leading index number.
 * \param data CCVF leading index, cell int data value (supplied).
 */
    // Set a cell dim value for a C++ CAR_CU_Mesh int CCVF class object
    // (self) - works for both the arbitrary leading index size and the 
    // default with the size of the leading index equal to that of the
    // problem geometry dimension.
    void set_mesh_ccvf_i_cell_dim_(long & mesh_ind, long & self, 
				   long & cell_ind, long & dim_ind, 
				   long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<int> CCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_ccvf_i_cell_dim_!");
	Insist(idim > 0 && idim <= CCVF_int.get_size(), 
	       "Invalid dim number passed to set_mesh_ccvf_i_cell_dim_!");

	CCVF_int(idim, icell) = static_cast<int>(data);

    }

//===========================================================================//
// double CCVF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double cell-centered vector field (CCVF) object 
 *        that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 * \param data CCVF double data values (returned)
 * \param data_size Size of the CCVF returned data vector (must equal the 
 *                  number of cells in the mesh times the leading index size).
 */
    // Return an entire C++ CAR_CU_Mesh double CCVF class object (self) - 
    // works for both the arbitrary leading index size and the default with 
    // the size of the leading index equal to that of the problem geometry 
    // dimension.
    void get_mesh_ccvf_d_(long & mesh_index, long & self, double & data, 
			  long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCVF<double> CCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size ==  mesh->get_num_cells() * CCVF_double.get_size(), 
	       "Invalid data size passed to get_mesh_ccvf_d_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int dim = 1; dim <= CCVF_double.get_size(); dim++)
	    {
	        * data_array = CCVF_double(dim, cell);
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the leading index values for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class double
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 * \param cell_ind The cell number.
 * \param data CCVF double cell data values (returned).
 * \param data_size Size of the CCVF returned data vector (must equal the 
 *                  size of the CCVF leading index).
 */
    // Return all of the dim values for a cell in a C++ CAR_CU_Mesh double 
    // CCVF class object (self) - works for both the arbitrary leading index 
    // size and the default with the size of the leading index equal to that 
    // of the problem geometry dimension.
    void get_mesh_ccvf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<double> CCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_ccvf_d_cell_!");
	Insist(idata_size == CCVF_double.get_size(), 
	       "Invalid data size passed to get_mesh_ccvf_d_cell_!");

	for (int idim = 1; idim <= CCVF_double.get_size(); idim++)
	{
	    * data_array = CCVF_double(idim, icell);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified leading index value for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class double 
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 * \param cell_ind The cell number.
 * \param dim_ind The leading index number.
 * \param data CCVF double leading index, cell data value (returned).
 */
    // Return a cell dim value for a C++ CAR_CU_Mesh double CCVF class
    // object (self) - works for both the arbitrary leading index size and 
    // the default with the size of the leading index equal to that of the 
    // problem geometry dimension.
    void get_mesh_ccvf_d_cell_dim_(long & mesh_ind, long & self, 
				   long & cell_ind, long & dim_ind,
				   double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<double> CCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_ccvf_d_cell_dim_!");
	Insist(idim > 0 && idim <= CCVF_double.get_size(), 
	       "Invalid dim number passed to get_mesh_ccvf_d_cell_dim_!");

	data = CCVF_double(idim, icell);

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double cell-centered vector field (CCVF) object
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 * \param data CCVF double data values (supplied).
 * \param data_size Size of the CCVF supplied data vector (must equal the 
 *                  number of cells in the mesh times the leading index size).
 */
    // Set an entire C++ CAR_CU_Mesh double CCVF class object (self) (can 
    // also be done at initialization using the constructor) - works for 
    // both the arbitrary leading index size and the default with the size 
    // of the leading index equal to that of the problem geometry dimension.
    void set_mesh_ccvf_d_(long & mesh_index, long & self, double & data, 
			  long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::CCVF<double> CCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size ==  mesh->get_num_cells() * CCVF_double.get_size(), 
	       "Invalid data size passed to set_mesh_ccvf_d_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int dim = 1; dim <= CCVF_double.get_size(); dim++)
	    {
	        CCVF_double(dim, cell) = * data_array;
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that sets all of the leading index values for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class double
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 * \param cell_ind The cell number.
 * \param data CCVF double cell data values (supplied).
 * \param data_size Size of the CCVF supplied data vector (must equal the 
 *                  size of the CCVF leading index).
 */
    // Set all of the dim values for a cell in a C++ CAR_CU_Mesh double 
    // CCVF class object (self) - works for both the arbitrary leading 
    // index size and the default with the size of the leading index equal 
    // to that of the problem geometry dimension.
    void set_mesh_ccvf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<double> CCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_ccvf_d_cell_!");
	Insist(idata_size == CCVF_double.get_size(), 
	       "Invalid data size passed to set_mesh_ccvf_d_cell_!");

	for (int idim = 1; idim <= CCVF_double.get_size(); idim++)
	{
	    CCVF_double(idim, icell) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified leading index value for the
 *        specified cell in the CAR_CU_Mesh nested mesh field class double
 *        cell-centered vector field (CCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double CCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double CCVF object.
 * \param cell_ind The cell number.
 * \param dim_ind The leading index number.
 * \param data CCVF double leading index, cell data value (supplied).
 */
    // Set a cell dim value for a C++ CAR_CU_Mesh double CCVF class object
    // (self) - works for both the arbitrary leading index size and the 
    // default with the size of the leading index equal to that of the 
    // problem geometry dimension.
    void set_mesh_ccvf_d_cell_dim_(long & mesh_ind, long & self, 
				   long & cell_ind, long & dim_ind,
				   double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and CCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::CCVF<double> CCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::CCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_ccvf_d_cell_dim_!");
	Insist(idim > 0 && idim <= CCVF_double.get_size(), 
	       "Invalid dim number passed to set_mesh_ccvf_d_cell_dim_!");

	CCVF_double(idim, icell) = data;

    }

//===========================================================================//
// int FCSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int face-centered scalar field (FCSF) object that
 *        is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCSF object.
 * \param data FCSF int data values (returned).
 * \param data_size Size of the FCSF returned data vector (must equal the 
 *                  number of cell faces in the mesh).
 */
    // Return an entire C++ CAR_CU_Mesh int FCSF class object (self).
    void get_mesh_fcsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<int> FCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes(), 
	       "Invalid data size passed to get_mesh_fcsf_i_!");

	for (int face = 1; face <= mesh->get_num_face_nodes(); face++)
	{
	    * data_array = static_cast<long>(FCSF_int(face));
	    ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns all of the face values for the specified
 *        cell in the CAR_CU_Mesh nested mesh field class int face-centered 
 *        scalar field (FCSF) object that is referenced by the specified 
 *        opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCSF object.
 * \param cell_ind The cell number. 
 * \param data FCSF int data values (returned).
 * \param data_size Size of the FCSF returned data vector (must equal the 
 *                  number of faces per cell).
 */
    // Return all of the face values for a cell in a C++ CAR_CU_Mesh integer
    // FCSF class object (self).
    void get_mesh_fcsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> FCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(),
	       "Invalid data size passed to get_mesh_fcsf_i_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = static_cast<long>(FCSF_int(icell, iface));
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified cell face value from the 
 *        CAR_CU_Mesh nested mesh field class int face-centered scalar field
 *        (FCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCSF int cell face data value (returned).
 */
    // Return a cell face value from a C++ CAR_CU_Mesh int FCSF class 
    // object (self).
    void get_mesh_fcsf_i_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind, 
				    long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> FCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_fcsf_i_cell_face_!");

	data = static_cast<long>(FCSF_int(icell, iface));

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int face-centered scalar field (FCSF) object that
 *        is referenced by the specified opaque pointers. This can also be 
 *        done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCSF object.
 * \param data FCSF int data values (supplied).
 * \param data_size Size of the FCSF supplied data vector (must equal the 
 *                  number of cell faces in the mesh).
 */
    // Set an entire C++ CAR_CU_Mesh int FCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<int> FCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes(), 
	       "Invalid data size passed to set_mesh_fcsf_i_!");

	for (int face = 1; face <= mesh->get_num_face_nodes(); face++)
	{
	    FCSF_int(face) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets all of the face values for the specified 
 *        cell in the CAR_CU_Mesh nested mesh field class int face-centered
 *        scalar field (FCSF) object that is referenced by the specified
 *        opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCSF object.
 * \param cell_ind The cell number. 
 * \param data FCSF cell face int data values (supplied).
 * \param data_size Size of the FCSF supplied data vector (must equal the 
 *                  number of faces per cell).
 */
    // Set all of the face values for a cell in a C++ CAR_CU_Mesh int FCSF
    // class object (self).
    void set_mesh_fcsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> FCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data size passed to set_mesh_fcsf_i_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    FCSF_int(icell, iface) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell face value in the 
 *        CAR_CU_Mesh nested mesh field class int face-centered scalar field
 *        (FCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCSF cell face int data value (supplied).
 */
    // Set a cell face value for a C++ CAR_CU_Mesh int FCSF class object
    // (self).
    void set_mesh_fcsf_i_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
			            long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<int> FCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to set_mesh_fcsf_i_cell_face_!");

	FCSF_int(icell, iface) = static_cast<int>(data);

    }

//===========================================================================//
// double FCSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double face-centered scalar field (FCSF) object
 *        that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 * \param data FCSF double data values (returned).
 * \param data_size Size of the FCSF returned data vector (must equal the 
 *                  number of cell faces in the mesh).
 */
    // Return an entire C++ CAR_CU_Mesh double FCSF class object (self).
    void get_mesh_fcsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<double> FCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes(), 
	       "Invalid data size passed to get_mesh_fcsf_d_!");

	for (int face = 1; face <= mesh->get_num_face_nodes(); face++)
	{
	    * data_array = FCSF_double(face);
	    ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns all of the face values for the specified
 *        cell in the CAR_CU_Mesh nested mesh field class double face-centered
 *        scalar field (FCSF) object that is referenced by the specified 
 *        opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 * \param cell_ind The cell number. 
 * \param data FCSF double data values (returned).
 * \param data_size Size of the FCSF returned data vector (must equal the 
 *                  number of faces per cell).
 */
    // Return all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCSF class object (self).
    void get_mesh_fcsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> FCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data size passed to get_mesh_fcsf_d_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = FCSF_double(icell, iface);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified cell face value from the 
 *        CAR_CU_Mesh nested mesh field class double face-centered scalar
 *        field (FCSF) object that is referenced by the specified opaque
 *        pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCSF double cell face data value (returned).
 */
    // Return a cell face value for a C++ CAR_CU_Mesh double FCSF class object
    // (self).
    void get_mesh_fcsf_d_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
				    double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> FCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_fcsf_d_cell_face_!");

	data = FCSF_double(icell, iface);

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double face-centered scalar field (FCSF) object
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 * \param data FCSF double data values (supplied).
 * \param data_size Size of the FCSF supplied data vector (must equal the 
 *                  number of cell faces in the mesh).
 */
    // Set an entire C++ CAR_CU_Mesh double FCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCSF<double> FCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes(), 
	       "Invalid data size passed to set_mesh_fcsf_d_!");

	for (int face = 1; face <= mesh->get_num_face_nodes(); face++)
	{
	    FCSF_double(face) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets all of the face values for the specified 
 *        cell in the CAR_CU_Mesh nested mesh field class double face-centered
 *        scalar field (FCSF) object that is referenced by the specified
 *        opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 * \param cell_ind The cell number. 
 * \param data FCSF cell face double data values (supplied).
 * \param data_size Size of the FCSF supplied data vector (must equal the 
 *                  number of faces per cell).
 */
    // Set all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCSF class object (self).
    void set_mesh_fcsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> FCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data size passed to set_mesh_fcsf_d_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    FCSF_double(icell, iface) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell face value in the 
 *        CAR_CU_Mesh nested mesh field class double face-centered scalar
 *        field (FCSF) object that is referenced by the specified opaque 
 *        pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCSF cell face double data value (supplied).
 */
    // Set a cell face value for a C++ CAR_CU_Mesh double FCSF class object
    // (self).
    void set_mesh_fcsf_d_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
				    double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCSF<double> FCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to set_mesh_fcsf_d_cell_face_!");

	FCSF_double(icell, iface) = data;

    }

//===========================================================================//
// int FCDSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int face-centered discontinuous scalar field (FCDSF)
 *        object that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object.
 * \param data FCSF int data values (returned).
 * \param data_size Size of the FCSF returned data vector (must equal the 
 *                  number of discontinuous cell faces in the mesh).
 */
    // Return an entire C++ CAR_CU_Mesh int FCDSF class object (self).
    void get_mesh_fcdsf_i_(long & mesh_index, long & self, 
			   long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<int> FCDSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data size passed to get_mesh_fcdsf_i_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        * data_array = static_cast<long>(FCDSF_int(cell, face));
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the face values for the specified
 *        cell in the CAR_CU_Mesh nested mesh field class int face-centered 
 *        discontinuous scalar field (FCDSF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object.
 * \param cell_ind The cell number. 
 * \param data FCDSF int data values (returned).
 * \param data_size Size of the FCDSF returned data vector (must equal the 
 *                  number of discontinuous faces per cell).
 */
     // Return all of the face values for a cell in a C++ CAR_CU_Mesh integer
    // FCDSF class object (self).
    void get_mesh_fcdsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			        long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> FCDSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcdsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(),
	       "Invalid data size passed to get_mesh_fcdsf_i_cell_!");

	vector<int> data_set = FCDSF_int(icell);

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = static_cast<long>(data_set[iface]);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified cell face value from the 
 *        CAR_CU_Mesh nested mesh field class int face-centered discontinuous
 *        scalar field (FCDSF) object that is referenced by the specified 
 *        opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCDSF int cell face data value (returned).
 */
    // Return a cell face value from a C++ CAR_CU_Mesh int FCDSF class 
    // object (self).
    void get_mesh_fcdsf_i_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind, 
				     long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> FCDSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcdsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_fcdsf_i_cell_face_!");

	data = static_cast<long>(FCDSF_int(icell, iface));

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int face-centered discontinuous scalar field 
 *        (FCDSF) object that is referenced by the specified opaque pointers.
 *        This can also be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object.
 * \param data FCDSF int data values (supplied).
 * \param data_size Size of the FCDSF supplied data vector (must equal the 
 *                  number of discontinuous cell faces in the mesh).
 */
    // Set an entire C++ CAR_CU_Mesh int FCDSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcdsf_i_(long & mesh_index, long & self, 
			   long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<int> FCDSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data size passed to set_mesh_fcdsf_i_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        FCDSF_int(cell, face) = static_cast<int>(* data_array);
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that sets all of the face values for the specified 
 *        cell in the CAR_CU_Mesh nested mesh field class int face-centered 
 *        discontinuous scalar field (FCDSF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object.
 * \param cell_ind The cell number. 
 * \param data FCDSF cell face int data values (supplied).
 * \param data_size Size of the FCDSF supplied data vector (must equal the 
 *                  number of discontinuous faces per cell).
 */
    // Set all of the face values for a cell in a C++ CAR_CU_Mesh int FCDSF
    // class object (self).
    void set_mesh_fcdsf_i_cell_(long & mesh_ind, long & self, long & cell_ind,
			        long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> FCDSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcdsf_i_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data size passed to set_mesh_fcdsf_i_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
            FCDSF_int(icell, iface) = static_cast<int>(* data_array);
     	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell face value in the 
 *        CAR_CU_Mesh nested mesh field class int face-centered discontinuous
 *        scalar field (FCDSF) object that is referenced by the specified
 *        opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCDSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCDSF cell face int data value (supplied).
 */
    // Set a cell face value for a C++ CAR_CU_Mesh int FCDSF class object
    // (self).
    void set_mesh_fcdsf_i_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind,
			             long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<int> FCDSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcdsf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to set_mesh_fcdsf_i_cell_face_!");

	FCDSF_int(icell, iface) = static_cast<int>(data);

    }

//===========================================================================//
// double FCDSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double face-centered discontinuous scalar field
 *        (FCDSF) object that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 * \param data FCSF double data values (returned).
 * \param data_size Size of the FCSF returned data vector (must equal the 
 *                  number of discontinuous cell faces in the mesh).
 */
    // Return an entire C++ CAR_CU_Mesh double FCDSF class object (self).
    void get_mesh_fcdsf_d_(long & mesh_index, long & self, 
			   double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<double> FCDSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size ==  mesh->get_num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data size passed to get_mesh_fcdsf_d_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        * data_array = FCDSF_double(cell, face);
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the face values for the specified
 *        cell in the CAR_CU_Mesh nested mesh field class double face-centered 
 *        discontinuous scalar field (FCDSF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 * \param cell_ind The cell number. 
 * \param data FCDSF double data values (returned).
 * \param data_size Size of the FCDSF returned data vector (must equal the 
 *                  number of discontinuous faces per cell).
 */
    // Return all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCDSF class object (self).
    void get_mesh_fcdsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			        double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> FCDSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcdsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data size passed to get_mesh_fcdsf_d_cell_!");

	vector<double> data_set = FCDSF_double(icell);

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    * data_array = data_set[iface];
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified cell face value from the 
 *        CAR_CU_Mesh nested mesh field class double face-centered 
 *        discontinuous scalar field (FCDSF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCDSF double cell face data value (returned).
 */
    // Return a cell face value for a C++ CAR_CU_Mesh double FCDSF class
    // object (self).
    void get_mesh_fcdsf_d_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind,
				     double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> FCDSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcdsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_fcdsf_d_cell_face_!");

	data = FCDSF_double(icell, iface);

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double face-centered discontinuous scalar field
 *        (FCDSF) object that is referenced by the specified opaque pointers.
 *        This can also be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 * \param data FCDSF double data values (supplied).
 * \param data_size Size of the FCDSF supplied data vector (must equal the 
 *                  number of discontinuous cell faces in the mesh).
 */
    // Set an entire C++ CAR_CU_Mesh double FCDSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcdsf_d_(long & mesh_index, long & self, 
			   double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCDSF<double> FCDSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size ==  mesh->get_num_cells() * 2 * mesh->get_ndim(), 
	       "Invalid data size passed to set_mesh_fcdsf_d_!");

	for (int cell = 1; cell <= mesh->get_num_cells(); cell++)
	{
	    for (int face = 1; face <= 2 * mesh->get_ndim(); face++)
	    {
	        FCDSF_double(cell, face) = * data_array;
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that sets all of the face values for the specified 
 *        cell in the CAR_CU_Mesh nested mesh field class double face-centered 
 *        discontinuous scalar field (FCDSF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 * \param cell_ind The cell number. 
 * \param data FCDSF cell face double data values (supplied).
 * \param data_size Size of the FCDSF supplied data vector (must equal the 
 *                  number of discontinuous faces per cell).
 */
    // Set all of the face values for a cell in a C++ CAR_CU_Mesh double 
    // FCDSF class object (self).
    void set_mesh_fcdsf_d_cell_(long & mesh_ind, long & self, long & cell_ind,
			        double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> FCDSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcdsf_d_cell_!");
	Insist(idata_size == 2 * mesh->get_ndim(), 
	       "Invalid data size passed to set_mesh_fcdsf_d_cell_!");

	for (int iface = 1; iface <= 2 * mesh->get_ndim(); iface++)
	{
	    FCDSF_double(icell, iface) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell face value in the 
 *        CAR_CU_Mesh nested mesh field class double face-centered 
 *        discontinuous scalar field (FCDSF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCDSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCDSF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCDSF cell face double data value (supplied).
 */
    // Set a cell face value for a C++ CAR_CU_Mesh double FCDSF class object
    // (self).
    void set_mesh_fcdsf_d_cell_face_(long & mesh_ind, long & self, 
				     long & cell_ind, long & face_ind,
				     double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCDSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCDSF<double> FCDSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::FCDSF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcdsf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to set_mesh_fcdsf_d_cell_face_!");

	FCDSF_double(icell, iface) = data;

    }

//===========================================================================//
// int FCVF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int face-centered vector field (FCVF) object
 *        that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCVF object.
 * \param data FCVF int data values (returned).
 * \param data_size Size of the FCVF returned data vector (must equal the 
 *                  number of unique cell faces in the mesh times the array
 *                  trailing index size).
 */
    // Return an entire C++ CAR_CU_Mesh int FCVF class object (self).
    void get_mesh_fcvf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCVF<int> fcvf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes() * fcvf_i.get_size(), 
	       "Invalid data size passed to get_mesh_fcvf_i_!");

	vector<vector<int> > data_set = fcvf_i();

	for (int node = 0; node < mesh->get_num_face_nodes(); node++)
	{
	    for (int dim = 0; dim < fcvf_i.get_size(); dim++)
	    {
		* data_array = static_cast<long>(data_set[node][dim]);
	        ++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the trailing index values for the
 *        specified cell face in the CAR_CU_Mesh nested mesh field class int 
 *        face-centered vector field (FCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCVF int data values (returned).
 * \param data_size Size of the FCVF returned data vector (must equal the 
 *                  size of the trailing index).
 */
    // Return all of the dim face values for a cell in a C++ CAR_CU_Mesh 
    // integer FCVF class object (self).
    void get_mesh_fcvf_i_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
				    long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<int> fcvf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcvf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_fcvf_i_cell_face_!");
	Insist(idata_size == fcvf_i.get_size(),
	       "Invalid data size passed to get_mesh_fcvf_i_cell_face_!");

	for (int idim = 1; idim <= fcvf_i.get_size(); idim++)
	{
	    * data_array = static_cast<long>(fcvf_i(icell, iface, idim));
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified trailing index cell face 
 *        value from the CAR_CU_Mesh nested mesh field class int face-centered 
 *        vector field (FCVF) object that is referenced by the specified 
 *        opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param dim_ind The trailing index number. 
 * \param data FCVF int cell face trailing index data value (returned).
 */
    // Return a dim cell face value from a C++ CAR_CU_Mesh int FCVF class 
    // object (self).
    void get_mesh_fcvf_i_cell_face_dim_(long & mesh_ind, long & self, 
					long & cell_ind, long & face_ind, 
					long & dim_ind, long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<int> fcvf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	    "Invalid cell number passed to get_mesh_fcvf_i_cell_face_dim_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	    "Invalid face number passed to get_mesh_fcvf_i_cell_face_dim_!");
	Insist(idim > 0 && idim <= fcvf_i.get_size(), 
	    "Invalid dimension passed to get_mesh_fcvf_i_cell_face_dim_!");

	data = static_cast<long>(fcvf_i(icell, iface, idim));

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int face-centered vector field (FCVF) object that
 *        is referenced by the specified opaque pointers. This can also be 
 *        done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCVF object.
 * \param data FCVF int data values (supplied).
 * \param data_size Size of the FCVF supplied data vector (must equal the 
 *                  number of unique cell faces in the mesh times the size 
 *                  of the trailing index).
 */
    // Set an entire C++ CAR_CU_Mesh int FCVF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcvf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCVF<int> fcvf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes() * fcvf_i.get_size(), 
	       "Invalid data size passed to set_mesh_fcvf_i_!");

	vector<vector<int> > data_set(mesh->get_num_face_nodes());

	for (int node = 0; node < mesh->get_num_face_nodes(); node++)
	{
	    data_set[node].resize(fcvf_i.get_size());
	    for (int dim = 0; dim < fcvf_i.get_size(); dim++)
	    {
	        data_set[node][dim] = static_cast<int>(* data_array);
	        ++data_array;
	    }
	}
	fcvf_i() = data_set;
    }

/*!
 * \brief Shadow object that sets all of the trailing index values for the 
 *        specified cell face in the CAR_CU_Mesh nested mesh field class int
 *        face-centered vector field (FCVF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCVF cell face int data values (supplied).
 * \param data_size Size of the FCVF supplied data vector (must equal the 
 *                  size of the vector trailing index).
 */
    // Set all of the dim face values for a cell in a C++ CAR_CU_Mesh int FCVF
    // class object (self).
    void set_mesh_fcvf_i_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind, 
				    long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<int> fcvf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcvf_i_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to set_mesh_fcvf_i_cell_face_!");
	Insist(idata_size == fcvf_i.get_size(),
	       "Invalid data size passed to set_mesh_fcvf_i_cell_face_!");

	for (int idim = 1; idim <= fcvf_i.get_size(); idim++)
	{
	    fcvf_i(icell, iface, idim) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell face trailing index value
 *        in the CAR_CU_Mesh nested mesh field class int face-centered vector
 *        field (FCVF) object that is referenced by the specified opaque
 *        pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param dim_ind The trailing index number. 
 * \param data FCVF cell face trailing index int data value (supplied).
 */
    // Set a dim cell face value for a C++ CAR_CU_Mesh int FCVF class object
    // (self).
    void set_mesh_fcvf_i_cell_face_dim_(long & mesh_ind, long & self, 
					long & cell_ind, long & face_ind,
					long & dim_ind, long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<int> fcvf_i = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<int> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	    "Invalid cell number passed to set_mesh_fcvf_i_cell_face_dim_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	    "Invalid face number passed to set_mesh_fcvf_i_cell_face_dim_!");
	Insist(idim > 0 && idim <= fcvf_i.get_size(), 
	    "Invalid dimension passed to get_mesh_fcvf_i_cell_face_dim_!");

	fcvf_i(icell, iface, idim) = static_cast<int>(data);

    }

//===========================================================================//
// double FCVF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double face-centered vector field (FCVF) object
 *        that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 * \param data FCVF double data values (returned).
 * \param data_size Size of the FCVF returned data vector (must equal the 
 *                  number of unique cell faces in the mesh times the array
 *                  trailing index size).
 */
    // Return an entire C++ CAR_CU_Mesh double FCVF class object (self).
    void get_mesh_fcvf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCVF<double> fcvf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes() * fcvf_d.get_size(), 
	       "Invalid data size passed to get_mesh_fcvf_d_!");

	vector<vector<double> > data_set = fcvf_d();

	for (int node = 0; node < mesh->get_num_face_nodes(); node++)
	{
	    for (int dim = 0; dim < fcvf_d.get_size(); dim++)
	    {
	        * data_array = data_set[node][dim];
	        ++data_array;
	    }
	}

    }

/*!
 * \brief Shadow object that returns all of the trailing index values for the
 *        specified cell face in the CAR_CU_Mesh nested mesh field class 
 *        double face-centered vector field (FCVF) object that is referenced
 *        by the specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCVF double data values (returned).
 * \param data_size Size of the FCVF returned data vector (must equal the 
 *                  size of the trailing index).
 */
    // Return all of the dim face values for a cell in a C++ CAR_CU_Mesh 
    // double FCVF class object (self).
    void get_mesh_fcvf_d_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind,
				    double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<double> fcvf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to get_mesh_fcvf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to get_mesh_fcvf_d_cell_face_!");
	Insist(idata_size == fcvf_d.get_size(),
	       "Invalid data size passed to get_mesh_fcvf_d_cell_face_!");

	for (int idim = 1; idim <= fcvf_d.get_size(); idim++)
	{
	    * data_array = fcvf_d(icell, iface, idim);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that returns the specified trailing index cell face 
 *        value from the CAR_CU_Mesh nested mesh field class double 
 *        face-centered vector field (FCVF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param dim_ind The trailing index number. 
 * \param data FCVF double cell face trailing index data value (returned).
 */
    // Return a dim cell face value from a C++ CAR_CU_Mesh double FCVF class 
    // object (self).
    void get_mesh_fcvf_d_cell_face_dim_(long & mesh_ind, long & self, 
					long & cell_ind, long & face_ind, 
					long & dim_ind,  double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<double> fcvf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	    "Invalid cell number passed to get_mesh_fcvf_d_cell_face_dim_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	    "Invalid face number passed to get_mesh_fcvf_d_cell_face_dim_!");
	Insist(idim > 0 && idim <= fcvf_d.get_size(), 
	    "Invalid dimension passed to get_mesh_fcvf_d_cell_face_dim_!");

	data = fcvf_d(icell, iface, idim);

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double face-centered vector field (FCVF) object
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 * \param data FCVF double data values (supplied).
 * \param data_size Size of the FCVF supplied data vector (must equal the 
 *                  number of unique cell faces in the mesh times the size 
 *                  of the trailing index).
 */
    // Set an entire C++ CAR_CU_Mesh double FCVF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_fcvf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::FCVF<double> fcvf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == mesh->get_num_face_nodes() * fcvf_d.get_size(), 
	       "Invalid data size passed to set_mesh_fcvf_d_!");

	vector<vector<double> > data_set(mesh->get_num_face_nodes());

	for (int node = 0; node < mesh->get_num_face_nodes(); node++)
	{
	    data_set[node].resize(fcvf_d.get_size());
	    for (int dim = 0; dim < fcvf_d.get_size(); dim++)
	    {
		data_set[node][dim] = * data_array;
		++data_array;
	    }
	}
	fcvf_d() = data_set;
    }

/*!
 * \brief Shadow object that sets all of the trailing index values for the 
 *        specified cell face in the CAR_CU_Mesh nested mesh field class double
 *        face-centered vector field (FCVF) object that is referenced by the
 *        specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param data FCVF cell face double data values (supplied).
 * \param data_size Size of the FCVF supplied data vector (must equal the 
 *                  size of the vector trailing index).
 */
    // Set all of the dim face values for a cell in a C++ CAR_CU_Mesh double 
    // FCVF class object (self).
    void set_mesh_fcvf_d_cell_face_(long & mesh_ind, long & self, 
				    long & cell_ind, long & face_ind, 
				    double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<double> fcvf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	       "Invalid cell number passed to set_mesh_fcvf_d_cell_face_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	       "Invalid face number passed to set_mesh_fcvf_d_cell_face_!");
	Insist(idata_size == fcvf_d.get_size(),
	       "Invalid data size passed to set_mesh_fcvf_d_cell_face_!");

	for (int idim = 1; idim <= fcvf_d.get_size(); idim++)
	{
	    fcvf_d(icell, iface, idim) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified cell face trailing index value
 *        in the CAR_CU_Mesh nested mesh field class double face-centered 
 *        vector field (FCVF) object that is referenced by the specified 
 *        opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double FCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double FCVF object.
 * \param cell_ind The cell number. 
 * \param face_ind The cell face number. 
 * \param dim_ind The trailing index number. 
 * \param data FCVF cell face trailing index double data value (supplied).
 */
    // Set a dim cell face value for a C++ CAR_CU_Mesh double FCVF class object
    // (self).
    void set_mesh_fcvf_d_cell_face_dim_(long & mesh_ind, long & self, 
					long & cell_ind, long & face_ind,
					long & dim_ind,  double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and FCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::FCVF<double> fcvf_d = 
	    * opaque_pointers<CAR_CU_Mesh::FCVF<double> >::item(self);
	// Cast the long variables to int
	int icell = static_cast<int>(cell_ind);
	int iface = static_cast<int>(face_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(icell > 0 && icell <= mesh->get_num_cells(), 
	    "Invalid cell number passed to set_mesh_fcvf_d_cell_face_dim_!");
	Insist(iface > 0 && iface <= 2 * mesh->get_ndim(), 
	    "Invalid face number passed to set_mesh_fcvf_d_cell_face_dim_!");
	Insist(idim  > 0 && idim <= fcvf_d.get_size(), 
	    "Invalid dimension passed to get_mesh_fcvf_d_cell_face_dim_!");

	fcvf_d(icell, iface, idim) = data;

    }

//===========================================================================//
// int NCSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int node-centered scalar field (NCSF) object that
 *        is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCSF object.
 * \param data NCSF int data values (returned).
 * \param data_size Size of the NCSF returned data vector (must equal 
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh).
 */
    // Return an entire C++ CAR_CU_Mesh int NCSF class object (self).
    void get_mesh_ncsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCSF<int> NCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == NCSF_int.get_size(), 
               "Invalid data size passed to get_mesh_ncsf_i_!");

	for (int node = 0; node < NCSF_int.get_size(); node++)
	{
	    * data_array = static_cast<long>(NCSF_int(node));
	    ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns the specified node value from the 
 *        CAR_CU_Mesh nested mesh field class int node-centered scalar field
 *       (NCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCSF object.
 * \param node_ind The node number. 
 * \param data NCSF int node data value (returned).
 */
    // Return a node value from a C++ CAR_CU_Mesh int NCSF class object
    // (self).
    void get_mesh_ncsf_i_node_(long & mesh_ind, long & self, long & node_ind,
			       long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCSF<int> NCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<int> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);

	Insist(inode > 0 && inode <= NCSF_int.get_size(), 
	       "Invalid node number passed to get_mesh_ncsf_i_node_!");

	data = static_cast<long>(NCSF_int(inode));

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int node-centered scalar field (NCSF) object 
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCSF object.
 * \param data NCSF int data values (supplied).
 * \param data_size Size of the NCSF supplied data vector (must equal 
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh).
 */
    // Set an entire C++ CAR_CU_Mesh int NCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ncsf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCSF<int> NCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == NCSF_int.get_size(), 
	       "Invalid data size passed to set_mesh_ncsf_i_!");

	for (int node = 0; node < NCSF_int.get_size(); node++)
	{
	    NCSF_int(node) = static_cast<int>(* data_array);
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified node value for the 
 *        CAR_CU_Mesh nested mesh field class int node-centered scalar field
 *        (NCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCSF object.
 * \param node_ind The node number. 
 * \param data NCSF node int data value (supplied).
 */
    // Set a node value for a C++ CAR_CU_Mesh int NCSF class object
    // (self).
    void set_mesh_ncsf_i_node_(long & mesh_ind, long & self, long & node_ind,
			       long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCSF<int> NCSF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<int> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);

	Insist(inode > 0 && inode <= NCSF_int.get_size(), 
		   "Invalid node number passed to set_mesh_ncsf_i_node_!");

	NCSF_int(inode) = static_cast<int>(data);

    }

//===========================================================================//
// double NCSF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double node-centered scalar field (NCSF) object that
 *        is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCSF object.
 * \param data NCSF double data values (returned).
 * \param data_size Size of the NCSF returned data vector (must equal
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh).
 */
    // Return an entire C++ CAR_CU_Mesh double NCSF class object (self).
    void get_mesh_ncsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCSF<double> NCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == NCSF_double.get_size(), 
	       "Invalid data size passed to get_mesh_ncsf_d_!");

	for (int node = 0; node < NCSF_double.get_size(); node++)
	{
	    * data_array = NCSF_double(node);
	    ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns the specified node value from the 
 *        CAR_CU_Mesh nested mesh field class double node-centered scalar field
 *       (NCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCSF object.
 * \param node_ind The node number. 
 * \param data NCSF double node data value (returned).
 */
    // Return a node value from a C++ CAR_CU_Mesh double NCSF class object
    // (self).
    void get_mesh_ncsf_d_node_(long & mesh_ind, long & self, long & node_ind,
			       double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCSF<double> NCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<double> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);

	Insist(inode > 0 && inode <= NCSF_double.get_size(), 
	       "Invalid node number passed to get_mesh_ncsf_d_node_!");

	data = NCSF_double(inode);

    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double node-centered scalar field (NCSF) object 
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCSF object.
 * \param data NCSF double data values (supplied).
 * \param data_size Size of the NCSF supplied data vector (must equal 
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh).
 */
    // Set an entire C++ CAR_CU_Mesh double NCSF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ncsf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCSF<double> NCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == NCSF_double.get_size(), 
	       "Invalid data size passed to set_mesh_ncsf_d_!");

	for (int node = 0; node < NCSF_double.get_size(); node++)
	{
	    NCSF_double(node) = * data_array;
	    ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified node value for the 
 *        CAR_CU_Mesh nested mesh field class double node-centered scalar field
 *        (NCSF) object that is referenced by the specified opaque pointers. 
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double NCSF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCSF object.
 * \param node_ind The node number. 
 * \param data NCSF node double data value (supplied).
 */
    // Set a node value for a C++ CAR_CU_Mesh double NCSF class object
    // (self).
    void set_mesh_ncsf_d_node_(long & mesh_ind, long & self, long & node_ind,
			       double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCSF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCSF<double> NCSF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCSF<double> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);

	Insist(inode > 0 && inode <= NCSF_double.get_size(), 
		   "Invalid node number passed to set_mesh_ncsf_d_node_!");

	NCSF_double(inode) = data;

    }

//===========================================================================//
// int NCVF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class int node-centered vector field (NCVF) object that
 *        is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this int NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCVF object.
 * \param data NCVF int data values (returned).
 * \param data_size Size of the NCVF returned data vector (must equal 
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh times the size of the trailing index).
 */
    // Return an entire C++ CAR_CU_Mesh int NCVF class object (self).
    void get_mesh_ncvf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCVF<int> NCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == NCVF_int.get_size_1() * NCVF_int.get_size_2(), 
               "Invalid data size passed to get_mesh_ncvf_i_!");

	for (int node = 1; node <= NCVF_int.get_size_1(); node++)
	{
	    for (int dim = 1; dim <= NCVF_int.get_size_2(); dim++)
	    {
	        * data_array = static_cast<long>(NCVF_int(node, dim));
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the trailing index values for the
 *        specified node in the CAR_CU_Mesh nested mesh field class int
 *        node-centered vector field (NCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCVF object.
 * \param node_ind The node number.
 * \param data NCVF int node data values (returned).
 * \param data_size Size of the NCVF returned data vector (must equal the 
 *                  size of the NCVF trailing index).
 */
    // Return all of the dim values for a node from a C++ CAR_CU_Mesh int NCVF
    // class object (self).
    void get_mesh_ncvf_i_node_(long & mesh_ind, long & self, long & node_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<int> NCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(inode > 0 && inode <= NCVF_int.get_size_1(), 
	       "Invalid node number passed to get_mesh_ncvf_i_node_!");
	Insist(idata_size == NCVF_int.get_size_2(), 
               "Invalid data size passed to get_mesh_ncvf_i_!");

	vector<int> data_set = NCVF_int(inode);

	for (int dim = 1; dim <= NCVF_int.get_size_2(); dim++)
	{
            * data_array = static_cast<long>(data_set[dim]);
            ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns the specified trailing index value for the
 *        specified node in the CAR_CU_Mesh nested mesh field class int 
 *        node-centered vector field (NCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCVF object.
 * \param node_ind The node number.
 * \param dim_ind The trailing index number.
 * \param data NCVF int trailing index, node data value (returned).
 */
    // Return a dim node value for a C++ CAR_CU_Mesh int NCVF class object
    // (self).
    void get_mesh_ncvf_i_node_dim_(long & mesh_ind, long & self, 
				   long & node_ind, long & dim_ind,
				   long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<int> NCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(inode > 0 && inode <= NCVF_int.get_size_1(), 
	       "Invalid node number passed to set_mesh_ncvf_i_node_dim!");
	Insist(idim> 0 && idim <=NCVF_int.get_size_2(), 
               "Invalid vector index passed to set_mesh_ncvf_i_node_dim_!");

	data = static_cast<long>(NCVF_int(inode, idim));
    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class int node-centered vector field (NCVF) object that
 *        is referenced by the specified opaque pointers. This can also be
 *        done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that
 *                   contains this int NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCVF object.
 * \param data NCVF int data values (supplied).
 * \param data_size Size of the NCVF supplied data vector  (must equal 
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh times the size of the trailing index)..
 */
    // Set an entire C++ CAR_CU_Mesh int NCVF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ncvf_i_(long & mesh_index, long & self, 
			  long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCVF<int> NCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(idata_size == NCVF_int.get_size_1() * NCVF_int.get_size_2(), 
	       "Invalid data size passed to set_mesh_ncvf_i_!");

	for (int node = 1; node <= NCVF_int.get_size_1(); node++)
	{
	    for (int dim = 1; dim <= NCVF_int.get_size_2(); dim++)
	    {
	        NCVF_int(node, dim) = static_cast<int>(* data_array);
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that sets all of the trailing index values for the
 *        specified node in the CAR_CU_Mesh nested mesh field class int 
 *        node-centered vector field (NCVF) object that is referenced by 
 *        the specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCVF object.
 * \param node_ind The node number.
 * \param data NCVF int node data values (supplied).
 * \param data_size Size of the NCVF supplied data vector (must equal the 
 *                  size of the NCVF trailing index).
 */
    // Set all of a nodes values for a C++ CAR_CU_Mesh int NCVF class object
    // (self).
    void set_mesh_ncvf_i_node_(long & mesh_ind, long & self, long & node_ind,
			       long & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<int> NCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idata_size = static_cast<int>(data_size);
	long * data_array = & data;

	Insist(inode > 0 && inode <= NCVF_int.get_size_1(), 
	       "Invalid node number passed to set_mesh_ncvf_i_node_!");
	Insist(idata_size == NCVF_int.get_size_2(), 
               "Invalid data size passed to set_mesh_ncvf_i_node_!");

	vector<int> data_set = NCVF_int(inode);

	for (int dim = 1; dim <= NCVF_int.get_size_2(); dim++)
	{
            data_set[dim] = static_cast<int>(* data_array);
            ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified trailing index value for the
 *        specified node in the CAR_CU_Mesh nested mesh field class int 
 *        node-centered vector field (NCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this int NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             int NCVF object.
 * \param node_ind The node number.
 * \param dim_ind The trailing index number.
 * \param data NCVF int trailing index, node data value (supplied).
 */
    // Set a dim node value for a C++ CAR_CU_Mesh int NCVF class object
    // (self).
    void set_mesh_ncvf_i_node_dim_(long & mesh_ind, long & self, 
				   long & node_ind, long & dim_ind,
				   long & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<int> NCVF_int = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<int> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(inode > 0 && inode <= NCVF_int.get_size_1(), 
	       "Invalid node number passed to set_mesh_ncvf_i_node_dim_!");
	Insist(idim> 0 && idim <= NCVF_int.get_size_2(), 
               "Invalid vector index passed to set_mesh_ncvf_i_node_dim_!");

	NCVF_int(inode, idim) = static_cast<int>(data);
    }

//===========================================================================//
// double NCVF class objects
//===========================================================================//
/*!
 * \brief Shadow object that returns the entire specified CAR_CU_Mesh nested
 *        mesh field class double node-centered vector field (NCVF) object
 *        that is referenced by the specified opaque pointers.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCVF object.
 * \param data NCVF double data values (returned).
 * \param data_size Size of the NCVF returned data vector (must equal
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh times the size of the trailing index).
 */
    // Return an entire C++ CAR_CU_Mesh double NCVF class object (self).
    void get_mesh_ncvf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCVF<double> NCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == NCVF_double.get_size_1() * NCVF_double.get_size_2(), 
               "Invalid data size passed to get_mesh_ncvf_d_!");

	for (int node = 1; node <= NCVF_double.get_size_1(); node++)
	{
	    for (int dim = 1; dim <= NCVF_double.get_size_2(); dim++)
	    {
	        * data_array = NCVF_double(node, dim);
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that returns all of the trailing index values for the
 *        specified node in the CAR_CU_Mesh nested mesh field class double
 *        node-centered vector field (NCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCVF object.
 * \param node_ind The node number.
 * \param data NCVF double node data values (returned).
 * \param data_size Size of the NCVF returned data vector (must equal the 
 *                  size of the NCVF trailing index).
 */
    // Return all of the dim values for a node from a C++ CAR_CU_Mesh double 
    // NCVF class object (self).
    void get_mesh_ncvf_d_node_(long & mesh_ind, long & self, long & node_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<double> NCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(inode > 0 && inode <= NCVF_double.get_size_1(), 
	       "Invalid node number passed to get_mesh_ncvf_d_node_!");
	Insist(idata_size == NCVF_double.get_size_2(), 
               "Invalid data size passed to get_mesh_ncvf_d_!");

	vector<double> data_set = NCVF_double(inode);

	for (int dim = 1; dim <= NCVF_double.get_size_2(); dim++)
	{
            * data_array = data_set[dim];
            ++data_array;
	}

    }

/*!
 * \brief Shadow object that returns the specified trailing index value for the
 *        specified node in the CAR_CU_Mesh nested mesh field class double 
 *        node-centered vector field (NCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCVF object.
 * \param node_ind The node number.
 * \param dim_ind The trailing index number.
 * \param data NCVF double trailing index, node data value (returned).
 */
    // Return a dim node value for a C++ CAR_CU_Mesh double NCVF class object
    // (self).
    void get_mesh_ncvf_d_node_dim_(long & mesh_ind, long & self, 
				   long & node_ind, long & dim_ind,
				   double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<double> NCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(inode > 0 && inode <= NCVF_double.get_size_1(), 
	       "Invalid node number passed to set_mesh_ncvf_d_node_dim!");
	Insist(idim> 0 && idim <=NCVF_double.get_size_2(), 
               "Invalid vector index passed to set_mesh_ncvf_d_node_dim_!");

	data = NCVF_double(inode, idim);
    }

/*!
 * \brief Shadow object that sets the entire specified CAR_CU_Mesh nested
 *        mesh field class double node-centered vector field (NCVF) object
 *        that is referenced by the specified opaque pointers. This can also
 *        be done at initialization using the constructor.
 * \param mesh_index Opaque pointer to the CAR_CU_Mesh class object that 
 *                   contains this double NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCVF object.
 * \param data NCVF double data values (supplied).
 * \param data_size Size of the NCVF supplied data vector(must equal
 *                  either the number of cell-corner plus face-centered nodes 
 *                  in the mesh or the number of cell-corner nodes in the 
 *                  mesh times the size of the trailing index).
 */
    // Set an entire C++ CAR_CU_Mesh double NCVF class object (self) (can 
    // also be done at initialization using the constructor).
    void set_mesh_ncvf_d_(long & mesh_index, long & self, 
			  double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_index) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_index);
	CAR_CU_Mesh::NCVF<double> NCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);
	// Cast the long variables to int
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(idata_size == NCVF_double.get_size_1() * NCVF_double.get_size_2(), 
	       "Invalid data size passed to set_mesh_ncvf_d_!");

	for (int node = 1; node <= NCVF_double.get_size_1(); node++)
	{
	    for (int dim = 1; dim <= NCVF_double.get_size_2(); dim++)
	    {
	        NCVF_double(node, dim) = * data_array;
		++data_array;
	    }
	}
    }

/*!
 * \brief Shadow object that sets all of the trailing index values for the
 *        specified node in the CAR_CU_Mesh nested mesh field class double 
 *        node-centered vector field (NCVF) object that is referenced by 
 *        the specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCVF object.
 * \param node_ind The node number.
 * \param data NCVF double node data values (supplied).
 * \param data_size Size of the NCVF supplied data vector (must equal the 
 *                  size of the NCVF trailing index).
 */
    // Set a node value for a C++ CAR_CU_Mesh double NCVF class object
    // (self).
    void set_mesh_ncvf_d_node_(long & mesh_ind, long & self, long & node_ind,
			       double & data, long & data_size)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<double> NCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idata_size = static_cast<int>(data_size);
	double * data_array = & data;

	Insist(inode > 0 && inode <= NCVF_double.get_size_1(), 
	       "Invalid node number passed to set_mesh_ncvf_d_node_!");
	Insist(idata_size == NCVF_double.get_size_2(), 
               "Invalid data size passed to set_mesh_ncvf_d_node_!");

	vector<double> data_set = NCVF_double(inode);

	for (int dim = 1; dim <= NCVF_double.get_size_2(); dim++)
	{
            data_set[dim] = * data_array;
            ++data_array;
	}
    }

/*!
 * \brief Shadow object that sets the specified trailing index value for the
 *        specified node in the CAR_CU_Mesh nested mesh field class double 
 *        node-centered vector field (NCVF) object that is referenced by the
 *        specified opaque pointers.
 * \param mesh_ind Opaque pointer to the CAR_CU_Mesh class object that 
 *                 contains this double NCVF class object.
 * \param self Opaque pointer to the CAR_CU_Mesh nested mesh field class 
 *             double NCVF object.
 * \param node_ind The node number.
 * \param dim_ind The trailing index number.
 * \param data NCVF double trailing index, node data value (supplied).
 */
    // Set a dim node value for a C++ CAR_CU_Mesh double NCVF class object
    // (self).
    void set_mesh_ncvf_d_node_dim_(long & mesh_ind, long & self, 
				   long & node_ind, long & dim_ind,
				   double & data)
    {
	// Get the addresses of the CAR_CU_Mesh (mesh_ind) and NCVF (self) 
        // class objects.
	SP<CAR_CU_Mesh> mesh = 
	    opaque_pointers<CAR_CU_Mesh>::item(mesh_ind);
	CAR_CU_Mesh::NCVF<double> NCVF_double = 
	    * opaque_pointers<CAR_CU_Mesh::NCVF<double> >::item(self);
	// Cast the long variables to int
	int inode = static_cast<int>(node_ind);
	int idim = static_cast<int>(dim_ind);

	Insist(inode > 0 && inode <= NCVF_double.get_size_1(), 
	       "Invalid node number passed to set_mesh_ncvf_d_node_dim_!");
	Insist(idim> 0 && idim <= NCVF_double.get_size_2(), 
               "Invalid vector index passed to set_mesh_ncvf_d_node_dim_!");

	NCVF_double(inode, idim) = data;
    }

} // end extern "C"

} // end namespace rtt_amr

#endif                          // __amr_Shadow_Mesh_cc__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Shadow_Mesh.cc
//---------------------------------------------------------------------------//
