//----------------------------------*-C++-*----------------------------------//
// XYCoord_sys.hh
// Thomas M. Evans
// Fri Jan 30 16:52:13 1998
//---------------------------------------------------------------------------//
// @> XYCoord_sys derived class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_XYCoord_sys_hh__
#define __imctest_XYCoord_sys_hh__

//===========================================================================//
// class XYCoord_sys - 
//
// Purpose : XY geometry coordinate system functions, derived
//           class of Coord_sys
//
// revision history:
// -----------------
//  0) original
//  1)  2-18-98 : added getCoord() function for debugging
//  2)  3-12-98 : moved Calc and Set_omega functions to Coord_sys because
//                they are the same for XY and XYZ transport, added transform 
//                for 2D meshes
//  3)  3-16-98 : added 
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Coord_sys.hh"
#include <vector>
#include <cmath>
#include <string>

IMCSPACE

using std::vector;
using std::sqrt;
using std::string;
    
class XYCoord_sys : public Coord_sys
{
public:
   // default constructor to set dimension of XY coordinate
  // system, inline
    XYCoord_sys() : Coord_sys(2) {}

  // virtual functions
    virtual string Get_coord() const { string c = "xy"; return c; }
    virtual double Transform(double dist_bnd, 
			     const vector<double> &omega) const 
    {
      // do transform back to xy plane if we know the distance to boundary on 
      // the xy plane
	return dist_bnd / sqrt(1.0 - omega[2] * omega[2]);
    }

  // pure virtual functions
    virtual vector<double> Calc_normal(vector<double> &center,
				       vector<double> &extent, int face)
    {
      // initialize normal vector (always 3D)
	vector<double> normal(Get_sdim());
	
      // step through faces to get the normal
	if (face == 1)
	{
	  // across -x face
	    normal[0] = -1.0;
	    normal[1] = 0.0;
	}
	else if (face == 2)
	{
	  // across +x face
	    normal[0] = 1.0;
	    normal[1] = 0.0;
	}
	else if (face == 3)
	{
	  // across -y face
	    normal[0] = 0.0;
	    normal[1] = -1.0;
	}
	else if (face == 4)
	{
	  // across +y face
	    normal[0] = 0.0;
	    normal[1] = 1.0;
	}

      // in XY geometry the z-normal is always 0
	normal[2] = 0.0;

      // return normal
	return normal;
    }
};

CSPACE

#endif                          // __imctest_XYCoord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/XYCoord_sys.hh
//---------------------------------------------------------------------------//
