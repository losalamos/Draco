//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Pyramid_Mesh.hh
 * \author Jeffery Densmore
 * \date   Mon Oct  6 09:15:12 2003
 * \brief  Pyramid_Mesh header file.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 *
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Pyramid_Mesh_hh
#define rtt_mc_Pyramid_Mesh_hh

#include "AMR_Layout.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "rng/Sprng.hh"
#include <vector>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Pyramid_Mesh

 * \brief An XYZ Pyramid mesh constructed from an R mesh.
 *
 * I'll put a description in later

 */

// revision history:
// -----------------
// 0) (Mon Oct  6 09:15:12 2003) Jeffery Densmore: original
// 
//===========================================================================//

class Pyramid_Mesh 
{

  public:
    // Typedefs used throughout Pyramid_Mesh class.
    typedef rtt_dsxx::SP<Coord_sys>           SP_Coord;
    typedef std::vector<double>               sf_double;
    typedef std::vector<int>                  sf_int;
    typedef std::vector<std::vector<double>>  vf_double;
    typedef std::string                       std_string;
    typedef std::pair<sf_double, sf_double>   pair_sf_double;
  
  private:
    // Base class reference to a derived coordinate system class
    // (the Pyramid_Mesh is always three-dimensional, Cartesian)
    SP_coord coord;

    // Layout (cell connectivity of mesh
    AMR_Layout layout;
    
    // vector<vector> of x-extents of each cell
    vf_double cell_x_extents;

    // Pyramid angle data; precalculated
    double beta_degrees;
    double beta_radians;
    double tan_beta;
    double sin_beta;
    double cos_beta;

    // Total, processor-local volume
    double total_volume;

    // Private implementations
    
    // Function to calculate frequently used wedge data
    void calc_angle_data(const double);

    // Function to calculate and set the total, on-processor volume
    void calc_total_volume();
  
  public:
    // Constructor
    Pyramid_Mesh(SP_Coord, AMR_Layout &, vf_double &, double);

    // Member functions used by Pyramid_Mesh dependent classes
    
    // Return the number of cells.
    int num_cells() const { return layout.num_cells(); }

    // Return the x-dimension cell extents
    double get_low_x(int cell) const {return cell_x_extents[cell-l][0]; }
    double get_high_x(int cell) const {return cell_x_extents[cell-1][1]; }

    // Determine if a position is in a cell
    bool in_cell(int, const sf_double &) const;

    // Services required for graphics dumps.
    sf_int get_cell_types() const;
    vf_double get_point_coord() const;

    // Required services for transport and source
    int get_cell(const sf_double *) const;
    double get_db(const sf_double &, const sf_double &, int, int &) const;
    inline sf_double get_normal(int,int) const;
    inline double volume(int) const;
    vf_double get_vertices(int) const;
    vf_double get_vertices(int, int) const;    
    int get_bndface(std_string, int) const;
    sf_int get_surcells(std_string) const;
    bool check_defined surcells(const std_string, const sf_int &) const;
    pair_sf_double sample_random_walk_sphere(int, const sf_double &, 
					     double, rng_Sprng &) const;



};
//---------------------------------------------------------------------------//
// PYRAMID_MESH INLINE FUNCTIONS
//---------------------------------------------------------------------------//



//---------------------------------------------------------------------------//
// RZWEDGE_MESH GENERALIZED MT SERVICES REQUIRED BY IMC
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the outward normal for the particular face of a cell.
 *
 * \param cell cell -- not used -- each Pyramid_Mesh cell has the same 
 * normals.
 * \param face face of cell for which to return outward normal
 *
 * \return outward normal
 */
Pyramid_Mesh::sf_double Pyramid_Mesh::get_normal(int cell, int face) const
{
    Check (coord->get_dim() == 3);
    Check ((face >=1) && (face<=6));

    sf_double normal(coord->get_dim(),0.0);

    // low x face
    if (face==1)
	normal[0]=-1.0;

    // high x face
    else if (face==2)
	normal[0]=1.0;

    // low y face
    else if (face==3)
    {
	normal[0]= -sin_beta;
	normal[1]= -cos_beta;
    }

    // high y face
    else if (face ==4)
    {
	normal[0]=-sin_beta;
	normal[1]=cos_beta;
    }
 

    // low z face
    else if (face==5)
    {
	normal[0]= -sin_beta;
	normal[2]= -cos_beta;
    }

    // high z face
    else if (face ==6)
    {
	normal[0]=-sin_beta;
	normal[2]=cos_beta;
    }
    
    // return outward normal
    return normal;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Return the volume of a Pyramid_Mesh cell.
 *
 * \param cell Pyramid_Mesh cell
 */
double Pyramid_Mesh::volume(int cell) const
{
    Require (cell>0 && cell <=num_cells());
    
    double lox = get_low_x(cell);
    double hix = get_high_x(cell);

    double vol = (4./3.)*tan_beta*tan_beta*(hix*hix*hix-lox*lox*lox);

    Ensure (vol>0.0);

    return vol;

}


} // end namespace rtt_mc

#endif // rtt_mc_Pyramid_Mesh_hh

//---------------------------------------------------------------------------//
//              end of mc/Pyramid_Mesh.hh
//---------------------------------------------------------------------------//
