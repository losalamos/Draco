//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/RZWedge_Builder.hh
 * \author Todd Urbatsch (extension of Thomas M. Evans's OS_Builder)
 * \date   Mon Feb  9 16:16:07 1998
 * \brief  RZWedge_Builder class header file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_RZWedge_Builder_hh__
#define __mc_RZWedge_Builder_hh__

#include "Coord_sys.hh"
#include "AMR_Layout.hh"
#include "RZWedge_Mesh.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class RZWedge_Builder
 *
 * The RZWedge_Builder class builds an instance of rtt_mc::RZWedge_Mesh using 
 * the same simple MC-defined orthogonal mesh format that OS_Mesh uses.  In
 * fact, RZWedge_Builder is merely a copy and extension of OS_Builder.  In the
 * future, the common parts of these two builders will be extracted.  This
 * less-than-ideal approach is taken to expedite milagro's use of an
 * RZWedge_Mesh.  
 * 
 */
// revision history:
// -----------------
//  0) original
//
//===========================================================================//

class RZWedge_Builder
{
  public:
    // Useful typdefs to std:: namespace members.
    typedef rtt_dsxx::SP<RZWedge_Mesh>        SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>           SP_Coord_sys;
    typedef rtt_dsxx::SP<AMR_Layout>          SP_Layout;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<int> >    vf_int;
    typedef std::vector<double>               sf_double;
    typedef std::vector<std::vector<double> > vf_double;
    typedef std::vector<std::string>          sf_string;
    typedef std::string                       std_string;
    typedef std::ifstream                     std_ifstream;
    
  private:
    // Data from Parser needed to build mesh.
    
    // Input file.
    std_string mesh_file;

    // Coordinate system string.
    std_string coord_system;

    // Wedge angle for 3D wedge representation of an RZ mesh.
    double theta_degrees;

    // Number of fine cells per coarse cell.
    vf_int fine_cells;
    
    // Ratio for fine zoning within a coarse cell (next fine cell 
    // is "ratio" times larger/smaller.)  Either all default to 
    // unity or all must be entered.
    vf_double fine_ratio;

    // Recursive total number of fine_cells per coarse cell.
    vf_int accum_cells;

    // Coarse edges.
    vf_double coarse_edge;

    // Number of fine_cells along each dimension.
    vf_double fine_edge;

    // Boundary conditions.
    sf_string bnd_cond;

    // Zone map.
    sf_int zone;
    vf_int cell_zone;
  
    // Defined cell regions.
    vf_int regions;

    // Surface source positional information.
    sf_string ss_pos;
    sf_int    num_defined_surcells;
    vf_int    defined_surcells;

    // Pointer to built Mesh.
    SP_Mesh mesh;

    // Member functions for building RZWedge_Mesh

    // Parse the mesh input file.
    void parser();
    void parser2D(std_ifstream &);
    void parser3D(std_ifstream &);
    void source_parser(std_ifstream &);

    // Build Layout helper functions.
    SP_Layout build_RZWedge_Layout(const Coord_sys &);
    void assignRZWedge_Layout(AMR_Layout &);

    // Build Coord_sys helper functions.
    SP_Coord_sys build_Coord();

    // Build Mesh helper functions.
    SP_Mesh build_RZWedge_Mesh(SP_Coord_sys, AMR_Layout &);

    // Member functions for cell-zone mapping
    void zone_mapper();
    void cell_zoner(int, int);
    void cell_zoner(int, int, int);

    // Calculate defined surface cells.
    void calc_defined_surcells();

  public:
    // Constructor.
    template<class IT> explicit RZWedge_Builder(rtt_dsxx::SP<IT>);

    // Build Mesh function.
    SP_Mesh build_Mesh();

    // Map a zone centered field into a cell centered field.
    template<class T> 
    std::vector<T> zone_cell_mapper(const std::vector<T> &) const;

    // ACCESSORS
    
    // Get a copy of the built mesh.
    SP_Mesh get_Mesh() const { Require(mesh); return mesh; }

    // Get cell regions for graphics dumping.
    sf_int get_regions() const;
    int get_num_regions() const { return regions.size(); }

    // Get cell zone information.
    int get_num_zones() const { return cell_zone.size(); }
    sf_int get_cells_in_zone(int z) const { return cell_zone[z-1]; }

    // Get the defined surcells list and positions.
    vf_int get_defined_surcells() const;
    sf_string get_ss_pos() const { return ss_pos; }

    //! Return the number of cells.
    int num_cells() const { return zone.size(); }
};

//---------------------------------------------------------------------------//
// Templated functions for RZWedge_Builder
//---------------------------------------------------------------------------//
// Constructor.

template<class IT>
RZWedge_Builder::RZWedge_Builder(rtt_dsxx::SP<IT> interface)
    : mesh_file(),
      coord_system(), 
      theta_degrees(),
      fine_cells(),
      fine_ratio(),
      accum_cells(),
      coarse_edge(),
      fine_edge(),
      bnd_cond(),
      zone(),
      cell_zone(),
      regions(),
      ss_pos(),
      num_defined_surcells(),
      defined_surcells(),
      mesh()
{
    Require (interface);

    // get mesh input file name from interface
    mesh_file = interface->get_mesh_file();

    // parse the mesh input file
    parser();

    // do zone mapping
    zone_mapper();

    Ensure (!mesh);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Map zone centered data fields into cell centered data fields.
 *
 * This function takes a zone centered field of size [0:Nz-1] where Nz is the
 * number of zones in the problem, and it maps the data into a cell centered
 * field of size [0:Nc-1].  Here Nc is the number of mesh cells in the
 * problem.  The builder has intrinsically (after construction) the
 * cell-to-zone mappings that make this operation possible.  The field types
 * must be vectors of type T.
 *
 * \param zone_field zone centered field to be converted into cell centered
 * field. 
 */
template<class T>
std::vector<T> RZWedge_Builder::zone_cell_mapper(const std::vector<T> &zone_field)
    const 
{
    // we will use vector throughout this function
    using std::vector;

    Require (zone_field.size() == cell_zone.size());

    // make a cell-sized return vector
    vector<T> cell_field(zone.size());

    // assign cell values to cell_field based on zonal values
    for (int cell = 0; cell < cell_field.size(); cell++)
	cell_field[cell] = zone_field[zone[cell]-1];

    // return the mapped field
    return cell_field;
}

} // end namespace rtt_mc

#endif                          // __mc_RZWedge_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/RZWedge_Builder.hh
//---------------------------------------------------------------------------//
