//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/OS_Builder.hh
 * \author Thomas M. Evans
 * \date   Mon Feb  9 16:16:07 1998
 * \brief  OS_Builder class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_OS_Builder_hh__
#define __mc_OS_Builder_hh__

#include "Coord_sys.hh"
#include "Layout.hh"
#include "OS_Mesh.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class OS_Builder
 *
 * The OS_Builder class builds an instance of rtt_mc::OS_Mesh using a simple
 * MC-defined orthogonal mesh format.  Future incantations may use more
 * advanced mesh format readers in order to build OS meshes from more
 * advanced formats.
 * 
 * /sa The examples page (mc/test/tstOSMesh.cc) for examples how to build OS
 * meshes using OS_Builder. 
 */
// revision history:
// -----------------
//  0) original
//  1)  3-18-98 : added generalized mesh constructor which consists of 
//                calculating vertex-based arrays and sending them to the 
//                OS_Mesh constructor
//  2)   4-6-99 : made OS_Builder class templated on an interface type
//  3)  4-13-99 : moved into mc package
//===========================================================================//

class OS_Builder
{
  public:
    // Useful typdefs to std:: namespace members.
    typedef dsxx::SP<OS_Mesh>                 SP_Mesh;
    typedef dsxx::SP<Coord_sys>               SP_Coord_sys;
    typedef dsxx::SP<Layout>                  SP_Layout;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<int> >    vf_int;
    typedef std::vector<double>               sf_double;
    typedef std::vector<std::vector<double> > vf_double;
    typedef std::vector<std::string>          sf_string;
    typedef std::string                       std_string;
    
  private:
    // Data from Parser needed to build mesh.

    // Coordinate system string.
    std_string coord_system;

    // Number of fine_cells along each dimension.
    vf_double fine_edge;

    // Boundary conditions.
    sf_string bnd_cond;
  
    // Member functions for building OS_Mesh

    // Build Layout helper functions.
    SP_Layout build_Layout(const Coord_sys &);
    void assign2D(Layout &);
    void assign3D(Layout &);

    // Build Coord_sys helper functions.
    SP_Coord_sys build_Coord();

    // Build Mesh helper functions.
    SP_Mesh build_2DMesh(SP_Coord_sys, Layout &);
    SP_Mesh build_3DMesh(SP_Coord_sys, Layout &);

  public:
    // Constructor.
    template<class IT> explicit OS_Builder(dsxx::SP<IT>);

    // Build Mesh function..
    SP_Mesh build_Mesh();
};

//---------------------------------------------------------------------------//
// inline functions for OS_Builder
//---------------------------------------------------------------------------//

template<class IT>
OS_Builder::OS_Builder(dsxx::SP<IT> interface)
{
    Require (interface);

    // get data arrays from OS_Interface needed to build OS_Mesh
    coord_system = interface->get_coordinates();
    fine_edge    = interface->get_fine_edge();
    bnd_cond     = interface->get_boundaries();
}

} // end namespace rtt_mc

#endif                          // __mc_OS_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/OS_Builder.hh
//---------------------------------------------------------------------------//
