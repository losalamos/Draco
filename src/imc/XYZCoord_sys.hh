//----------------------------------*-C++-*----------------------------------//
// XYZCoord_sys.hh
// Thomas M. Evans
// Fri Jan 30 16:45:36 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __imctest_XYZCoord_sys_hh__
#define __imctest_XYZCoord_sys_hh__

//===========================================================================//
// class XYZCoord_sys - 
//
// Date created : 1-12-98
// Purpose      : XYZ geometry coordinate system functions, derived
//                class of Coord_sys
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "Names.hh"
#include "Coord_sys.hh"
#include <vector>

IMCSPACE

using std::vector;
    
class XYZCoord_sys : public Coord_sys
{
public:
  // default constructor to set dimension of XY coordinate
  // system, inline
    XYZCoord_sys() : Coord_sys(3) {}
  // virtual functions
    virtual void setOmega(vector<double> &, Random &) const;
    virtual void calcOmega(double, double, vector<double> &) const;
};

CSPACE

#endif                          // __imctest_XYZCoord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/XYZCoord_sys.hh
//---------------------------------------------------------------------------//
