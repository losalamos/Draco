//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphyramid_Builder.hh
 * \author Jeffery Densmore (stolen from RZWedge_Builder.hh)
 * \date   Mon Nov  10 7:39:00 2003
 * \brief  Sphyramid_Builder class header file
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Sphyramid_Builder_hh__
#define __mc_Sphyramid_Builder_hh__

#include "Coord_sys.hh"
#include "Layout.hh"
#include "Sphyramid_Mesh.hh"
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
 * \class Sphyramid_Builder
 * \brief
 *
 * This class builds an instance of rtt_mc::Sphyramid_Mesh.  It is basically
 * a copy of RZWedge_Builder, with member functions changed to reflect new
 * Sphyramid geometry.
 *
 */
// revision history:
// -----------------
// 0) (Mon Nov  10 7:39:00 2003) Jeffery Densmore: original
// 
//===========================================================================//

class Sphyramid_Builder 
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Sphyramid_Mesh>        SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>             SP_Coord_sys;
    typedef rtt_dsxx::SP<Layout>                SP_Layout;
    typedef std::vector<int>                    sf_int;
    typedef std::vector<std::vector<int> >      vf_int;
    typedef std::vector<double>                 sf_double;
    typedef std::vector<std::vector<double> >   vf_double;
    typedef std::vector<std::string>            sf_string;
    typedef std::string                         std_string;
    typedef std::ifstream                       std_ifstream;

  private:
    // Data from parser needed to build mesh.

    // Input file.
    std_string mesh_file;

    // Coordinate system string.
    std_string coord_system;

    // Angle of spherical cone in degrees (i.e. maximum polar angle)
    double alpha_degrees;
    
    // Number of fine cells per coarse cell.
    sf_int fine_cells;

    // Ratio for fine zoning within a coarse cell (next fine cell is "ratio"
    // times large/smaller). Either all default to unity or all must be
    // entered.
    sf_double fine_ratio;

    // Total number of fine_cells per coarse cell.
    sf_int accum_cells;

    // Coarse edges.
    sf_double coarse_edge;

    // Fine edges.
    sf_double fine_edge;

    // Boundary conditions.
    sf_string bnd_cond;

    // Zone map.
    sf_int zone;
    vf_int cell_zone;

    //Defined cell regions
    vf_int regions;

    //Surface source positional information
    sf_string ss_pos;
    sf_int num_defined_surcells;
    vf_int defined_surcells;

    // Pointer to built Mesh.
    SP_Mesh mesh;

    // Member functions for building Sphyramid_Mesh

    // Parse the mesh input file.
    void parser();
    void parser1D(std_ifstream &in);
    void source_parser(std_ifstream &in);

    // Build Layout helper functions.
    SP_Layout build_Sphyramid_Layout(const Coord_sys &coord) const;
    void assign_Sphyramid_Layout(Layout &layout) const;

    // Build Coord_sys helper functions.
    SP_Coord_sys build_Coord() const;

    // Build Mesh helper functions
    SP_Mesh build_Sphyramid_Mesh(SP_Coord_sys coord, Layout &layout);

    // Member functions for cell-zone mapping
    void zone_mapper();
    void cell_zoner(int iz);

    // calculate defined surface cells.
    void calc_defined_surcells();

    // hide copy constructor and assignment operators
    Sphyramid_Builder(const Sphyramid_Builder &);
    Sphyramid_Builder & operator=(const Sphyramid_Builder &);

  public:
    //constructor
    template<class IT> explicit Sphyramid_Builder(rtt_dsxx::SP<IT> interface);

    //Build Mesh function.
    SP_Mesh build_Mesh();

    // Map a zone centered field into a cell centered field.
    template<class T>
    std::vector<T> zone_cell_mapper(const std::vector<T> &zone_field) const;

    // ACCESSORS

    // Get a copy of the built mesh.
    SP_Mesh get_Mesh() const { Require (mesh); return mesh; }

    // Get cell regions for graphics dumping.
    sf_int get_regions() const;
    int get_num_regions() const {return regions.size();}

    // Get cell zone information
    int get_num_zones() const { return cell_zone.size(); }
    sf_int get_cells_in_zone(int z) const { return cell_zone[z-1]; }

    // Get the defined surcells list of positions.
    vf_int get_defined_surcells() const;
    sf_string get_ss_pos() const { return ss_pos; }

    // Return the number of cells.
    int num_cells() const { return zone.size(); }
};

//---------------------------------------------------------------------------//
// Templated functions for Sphyramid_Builder
//---------------------------------------------------------------------------//
// Constructor.

template<class IT>
Sphyramid_Builder::Sphyramid_Builder(rtt_dsxx::SP<IT> interface)
    :mesh_file(),
     coord_system(),
     alpha_degrees(),
     fine_cells(),
     fine_ratio(),
     accum_cells(),
     coarse_edge(),
     fine_edge(),
     bnd_cond(2), // only 2 boundary conditions required in Sphyramid mesh
     zone(),
     cell_zone(),
     regions(),
     ss_pos(),
     num_defined_surcells(),
     defined_surcells(),
     mesh()
{
    Require (interface);

    // get mesh_input file name from interface
    this->mesh_file = interface->get_mesh_file();

    // parse the mesh input file
    parser();

    // do zone mapping
    zone_mapper();

    Ensure (!this->mesh);
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Map zone centered data fields into cell centered data fiels.
 *
 * This function takes a zone centered field of size[0:Nz-1] where Nz is the
 * number of zones in the problem, and it maps the data into a cell centered
 * field of size [0:Nc-1].  Here Nc is the number of mesh cells in the
 * problem.  The builder has intrinsically (after contstruction) the
 * cell-to-zone mappings that make this operation possible.  The field types
 * must be vectors of type T.
 * 
 * \param zone_field zone centered field to be converted into cell centered
 * field.
 *
 * \return cell-sized vector field
 */
template<class T>
std::vector<T> Sphyramid_Builder::zone_cell_mapper(const std::vector<T> 
						  &zone_field) const
{
    // we will use vector throughout this function
    using std::vector;

    Require (zone_field.size() == cell_zone.size());

    // make a cell-sized return vector
    vector<T> cell_field(zone.size());

    // assign cell values to cell_field based on zonal values
    for (int cell=0; cell < cell_field.size(); cell++)
	cell_field[cell] = zone_field[zone[cell]-1];

    // return the mapped field
    return cell_field;
}

} // end namespace rtt_mc

#endif // __mc_Sphyramid_Builder_hh__

//---------------------------------------------------------------------------//
//              end of mc/Sphyramid_Builder.hh
//---------------------------------------------------------------------------//
