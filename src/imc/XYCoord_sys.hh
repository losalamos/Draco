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
// 
//===========================================================================//

#include "imctest/Names.hh"
#include "imctest/Coord_sys.hh"
#include <vector>

IMCSPACE

using std::vector;
    
class XYCoord_sys : public Coord_sys
{
public:
  // default constructor to set dimension of XY coordinate
  // system, inline
    XYCoord_sys() : Coord_sys(2) {}
  // virtual functions
    virtual string Get_coord() const { string c = "xy"; return c; }
    virtual void Set_omega(vector<double> &, Random &) const;
    virtual void Calc_omega(double, double, vector<double> &) const;
};

CSPACE

#endif                          // __imctest_XYCoord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/XYCoord_sys.hh
//---------------------------------------------------------------------------//
