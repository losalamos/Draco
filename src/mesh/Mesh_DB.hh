//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/Mesh_DB.hh
 * \author Chris Gesh
 * \date   Mon Mar 13 08:09:08 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mesh_Mesh_DB_hh__
#define __mesh_Mesh_DB_hh__

#include "ds++/Mat.hh"
#include "nml/Group.hh"
#include "nml/Item.hh"

namespace rtt_mesh
{
 
//===========================================================================//
/*!
 * \class Mesh_DB
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Mesh_DB 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA
    
  public:
    
    void setup_namelist( NML_Group& g );

    // CREATORS

    // Public data to imitate NML

    int ncx;
    int ncy;
    int ncz;

    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;

    rtt_dsxx::Mat1<double> dx;
    rtt_dsxx::Mat1<double> dy;
    rtt_dsxx::Mat1<double> dz;

    Mesh_DB(int _ncx=10, int _ncy=10, int _ncz=10,
	    double xmin_in=0.0, double xmax_in=1.0, double ymin_in=0.0,
	    double ymax_in=1.0, double zmin_in=0.0, double zmax_in=1.0);

    Mesh_DB(const rtt_dsxx::Mat1<double> &dx_in,
	    const rtt_dsxx::Mat1<double> &dy_in,
    	    const rtt_dsxx::Mat1<double> &dz_in,
	    int ncx_in, int ncy_in, int ncz_in,
    	    double xmin_in=0.0, double xmax_in=1.0, double ymin_in=0.0,
    	    double ymax_in=1.0, double zmin_in=0.0, double zmax_in=1.0);

    //Mesh_DB(const Mesh_DB &rhs);
    //    ~Mesh_DB();

    // MANIPULATORS
    
    void resize();

    // Mesh_DB& operator=(const Mesh_DB &rhs);

    // ACCESSORS

    bool isValid() const;

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_mesh

#endif                          // __mesh_Mesh_DB_hh__

//---------------------------------------------------------------------------//
//                              end of mesh/Mesh_DB.hh
//---------------------------------------------------------------------------//
