//----------------------------------*-C++-*----------------------------------//
// XYZCoord_sys.hh
// Thomas M. Evans
// Fri Jan 30 16:45:36 1998
//---------------------------------------------------------------------------//
// @> XYZCoord_sys derived class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_XYZCoord_sys_hh__
#define __imctest_XYZCoord_sys_hh__

//===========================================================================//
// class XYZCoord_sys - 
//
// Purpose : XYZ geometry coordinate system functions, derived
//           class of Coord_sys
//
// revision history:
// -----------------
//  0) original
//  1)  2-18-98 : added getCoord() function for debugging
//  2)  3-12-98 : moved Calc and Set_omega functions to Coord_sys because
//                they are the same for XY and XYZ transport
//  3)  3-16-98 : reserve calc_normal function for later if need be
//  4)   5-5-98 : added sample_pos virtual function
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Coord_sys.hh"
#include "rng/Sprng.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>

IMCSPACE
    
using std::vector;
using std::string;

using RNG::Sprng;

class XYZCoord_sys : public Coord_sys
{
  // Begin_Doc xyzcoord_sys-int.tex
  // Begin_Verbatim 

public:
  // default constructor for 3D meshes
    XYZCoord_sys() : Coord_sys(3) {}

  // virtual functions
    virtual string get_Coord() const { string c = "xyz"; return c; }
    inline virtual vector<double> sample_pos(string, vector<double>,
					     vector<double>, Sprng &);
	     
  // End_Verbatim 
  // End_Doc 
};

//---------------------------------------------------------------------------//
// INLINE Functions
//---------------------------------------------------------------------------//
// sample the position in an XYZ cell

inline vector<double> 
XYZCoord_sys::sample_pos(string dist, vector<double> min, vector<double> max,
			 Sprng &random)
{
  // make return vector
    vector<double> r(3);

  // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);

    for (int d = 0; d < 3; d++)
    {
      // do uniform sampling
	if (dist == "uniform")
	    r[d] = (max[d] - min[d]) * random.ran() + min[d];
      // add others as needed
    }

  // return assigned array
    return r;
}

CSPACE

#endif                          // __imctest_XYZCoord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/XYZCoord_sys.hh
//---------------------------------------------------------------------------//
