//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/IMC_Test.cc
 * \author Thomas M. Evans
 * \date   Thu Aug  3 14:20:15 2000
 * \brief  Definitions from IMC_Test.hh.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "../Global.hh"
#include "mc/AMR_Layout.hh"
#include "mc/XYZCoord_sys.hh"
#include <cmath>

namespace rtt_imc_test
{

typedef std::vector<double>               sf_double;
typedef std::vector<std::vector<double> > vf_double;

//===========================================================================//
// BUILD A RZWEDGE MESH
//===========================================================================//

rtt_dsxx::SP<rtt_mc::RZWedge_Mesh> make_RZWedge_Mesh_AMR(double phi)
{
    using rtt_mc::global::pi;
    using rtt_mc::RZWedge_Mesh;
    using rtt_mc::Coord_sys;
    using rtt_mc::XYZCoord_sys;
    using rtt_mc::AMR_Layout;
    using rtt_dsxx::SP;
    using std::sqrt;
    using std::sin;
    using std::cos;

    double rphi = phi * pi / 180.0;
    
    // make coord
    SP<Coord_sys> coord(new XYZCoord_sys());

    // make layout (single level mesh)
    AMR_Layout layout(9, 6);
    {
	layout(1, 1, 1) = 1;
	layout(1, 2, 1) = 4;
	layout(1, 3, 1) = 1;
	layout(1, 4, 1) = 1;
	layout(1, 5, 1) = 1;
	layout(1, 6, 1) = 2;
    
	layout(2, 1, 1) = 2;
	layout(2, 2, 1) = 5;
	layout(2, 3, 1) = 2;
	layout(2, 4, 1) = 2;
	layout(2, 5, 1) = 1;
	layout(2, 6, 1) = 3;
    
	layout.set_size(3, 2, 2);
	layout(3, 1, 1) = 3;
	layout(3, 2, 1) = 6;
	layout(3, 2, 2) = 7;
	layout(3, 3, 1) = 3;
	layout(3, 4, 1) = 3;
	layout(3, 5, 1) = 2;
	layout(3, 6, 1) = 0;
    
	layout(4, 1, 1) = 1;
	layout(4, 2, 1) = 0;
	layout(4, 3, 1) = 4;
	layout(4, 4, 1) = 4;
	layout(4, 5, 1) = 4;
	layout(4, 6, 1) = 5;
    
	layout.set_size(5, 6, 2);
	layout(5, 1, 1) = 2;
	layout(5, 2, 1) = 0;
	layout(5, 3, 1) = 5;
	layout(5, 4, 1) = 5;
	layout(5, 5, 1) = 4;
	layout(5, 6, 1) = 6;
	layout(5, 6, 2) = 8;
    
	layout(6, 1, 1) = 3;
	layout(6, 2, 1) = 8;
	layout(6, 3, 1) = 6;
	layout(6, 4, 1) = 6;
	layout(6, 5, 1) = 5;
	layout(6, 6, 1) = 7;
    
	layout(7, 1, 1) = 3;
	layout(7, 2, 1) = 9;
	layout(7, 3, 1) = 7;
	layout(7, 4, 1) = 7;
	layout(7, 5, 1) = 6;
	layout(7, 6, 1) = 0;
    
	layout(8, 1, 1) = 6;
	layout(8, 2, 1) = 0;
	layout(8, 3, 1) = 8;
	layout(8, 4, 1) = 8;
	layout(8, 5, 1) = 5;
	layout(8, 6, 1) = 9;
    
	layout(9, 1, 1) = 7;
	layout(9, 2, 1) = 0;
	layout(9, 3, 1) = 9;
	layout(9, 4, 1) = 9;
	layout(9, 5, 1) = 8;
	layout(9, 6, 1) = 0;
    }

    // make x-z extents for each cell
    vf_double xz(9, sf_double(4));
    vf_double rz(9, sf_double(4));
    {
	double rconv = std::sqrt(rphi/std::sin(rphi))*std::cos(rphi/2);

	rz[0][0] = 0/rconv;
	rz[0][1] = 1/rconv;
	rz[0][2] = 0;
	rz[0][3] = 1;
	
	rz[1][0] = 0/rconv;
	rz[1][1] = 1/rconv;
	rz[1][2] = 1;
	rz[1][3] = 2;
	
	rz[2][0] = 0/rconv;
	rz[2][1] = 1/rconv;
	rz[2][2] = 2;
	rz[2][3] = 3;
	
	rz[3][0] = 1/rconv;
	rz[3][1] = 2/rconv;
	rz[3][2] = 0;
	rz[3][3] = 1;
	
	rz[4][0] = 1/rconv;
	rz[4][1] = 2/rconv;
	rz[4][2] = 1;
	rz[4][3] = 2;
	
	rz[5][0] = 1/rconv;
	rz[5][1] = 1.5/rconv;
	rz[5][2] = 2;
	rz[5][3] = 2.5;
	
	rz[6][0] = 1/rconv;
	rz[6][1] = 1.5/rconv;
	rz[6][2] = 2.5;
	rz[6][3] = 3;
	
	rz[7][0] = 1.5/rconv;
	rz[7][1] = 2/rconv;
	rz[7][2] = 2;
	rz[7][3] = 2.5;
	
	rz[8][0] = 1.5/rconv;
	rz[8][1] = 2/rconv;
	rz[8][2] = 2.5;
	rz[8][3] = 3;

	for (int i = 0; i < layout.num_cells(); i++)
	{
	    // convert rz faces
	    for (int j = 0; j < 2; j++)
	    {
		xz[i][j]   = sqrt(rphi/sin(rphi)) * rz[i][j] * cos(rphi/2.0); 
		xz[i][j+2] = rz[i][j+2];
	    }
	}
    }

    // make the mesh
    SP<RZWedge_Mesh> mesh(new RZWedge_Mesh(coord, layout, xz, phi));
    Check (mesh->full_Mesh());

    return mesh;
}

} // end namespace rtt_imc_test


//---------------------------------------------------------------------------//
//                              end of IMC_Test.cc
//---------------------------------------------------------------------------//
