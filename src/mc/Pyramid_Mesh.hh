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
    typedef std::vector<std::vector<double>>  vf_double;
  
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

    // Function to calculate and set the total, on-processor voluem
    void calc_total_volume();
  
  public:
    // Constructor
    Pyramid_Mesh(SP_Coord, AMR_Layout &, vf_double &, double);


  



};

} // end namespace rtt_mc

#endif // rtt_mc_Pyramid_Mesh_hh

//---------------------------------------------------------------------------//
//              end of mc/Pyramid_Mesh.hh
//---------------------------------------------------------------------------//
