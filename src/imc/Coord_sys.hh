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
// 0) original
// 
//===========================================================================//

#include "Names.hh"
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
public:
  // constructor for setting dimension of Coord_sys, inline
    Coord_sys(int dimension_) 
	:dimension(dimension_) {}
    int getDim() const { return dimension; } 
  // virtual functions utilized by each derived class
    virtual string getCoord() const = 0;
    virtual void setOmega(vector<double> &, Random &) const = 0;
    virtual void calcOmega(double, double, vector<double> &) const = 0;
};

CSPACE

#endif                          // __imctest_Coord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Coord_sys.hh
//---------------------------------------------------------------------------//
