//----------------------------------*-C++-*----------------------------------//
// Coord_sys.hh
// Thomas M. Evans
// Fri Jan 30 16:36:51 1998
//---------------------------------------------------------------------------//
// @> Coord_sys abstract base class header file
//---------------------------------------------------------------------------//

#ifndef __imctest_Coord_sys_hh__
#define __imctest_Coord_sys_hh__

//===========================================================================//
// class Coord_sys - 
//
// Purpose : abstract base class which defines the coordinate system
//           which the Mesh lives in
//
// revision history:
// -----------------
//  0) original
//  1)  3-12-98 : moved Calc and Set_omega functions into Coord_sys as
//                non-pure virtual functions because they are the same in
//                both XY and XYZ transport, added a transform for 2D meshes
// 
//===========================================================================//

#include "imctest/Names.hh"
#include <vector>
#include <string>

IMCSPACE

using std::vector;
using std::string;

class Random;

class Coord_sys
{
private:
  // dimension of system
    const int dimension;
    const int set_dimension;
public:
  // constructor for setting dimension of Coord_sys, inline
    Coord_sys(int dimension_) 
	:dimension(dimension_), set_dimension(3) {}

  // base class member functions

  // we have two dimensionalities, a "real" dimension for the geometry and a
  // "transport" dimension for MC transport which is inherently 3D
    int Get_dim() const { return dimension; } 
    int Get_sdim() const { return set_dimension; }


  // pure virtual functions
    virtual string Get_coord() const = 0;

  // virtual functions, these are only needed in some coordinate systems
    virtual double Transform(double dist_bnd, const vector<double> &omega) 
	const {	return dist_bnd; }

    virtual void Set_omega(vector<double> &, Random &) const;
    virtual void Calc_omega(double, double, vector<double> &) const;
};

CSPACE

#endif                          // __imctest_Coord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Coord_sys.hh
//---------------------------------------------------------------------------//
