//----------------------------------*-C++-*----------------------------------//
// Coord_sys.hh
// Thomas M. Evans
// Fri Jan 30 16:36:51 1998
//---------------------------------------------------------------------------//
// @> Coord_sys abstract base class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Coord_sys_hh__
#define __imc_Coord_sys_hh__

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
//  2)  3-16-98 : reserve calc_normal function for later if need be
//  3)  3-17-98 : because of a dumb-ass oversight on my part, we don't need
//                a transform for 2D XY, it has been removed
// 
//===========================================================================//

#include "imc/Names.hh"
#include "rng/Sprng.hh"
#include <vector>
#include <string>

IMCSPACE

using std::vector;
using std::string;
using RNG::Sprng;

class Coord_sys
{
private:
  // dimension of system
    const int dimension;
    const int set_dimension;

  // Begin_Doc coord_sys-int.tex 
  // Begin_Verbatim 

public:
  // constructor for setting dimension of Coord_sys, inline
    Coord_sys(int dimension_) 
	:dimension(dimension_), set_dimension(3) {}

  // virtual destructor to insure correct behavior down inheritance chain
    virtual ~Coord_sys() {}

  // base class member functions

  // we have two dimensionalities, a "real" dimension for the geometry and a
  // "transport" dimension for MC transport which is inherently 3D
    int get_dim() const { return dimension; } 
    int get_sdim() const { return set_dimension; }

  // pure virtual functions
    virtual string get_Coord() const = 0;
 
    virtual vector<double> 
    sample_pos(vector<double> &, vector<double> &, Sprng &) const = 0;

    virtual vector<double> 
    sample_pos(vector<double> &, vector<double> &, Sprng &, 
	       vector<double> &, double) const = 0;

    virtual 
    vector<double> sample_pos_on_face(vector<double> &, vector<double> &, 
				      int, Sprng &) const = 0;

  // virtual functions
    virtual vector<double> sample_dir(string, Sprng &) const;
    virtual void calc_omega(double, double, vector<double> &) const;


  // End_Verbatim 
  // End_Doc 
};

CSPACE

#endif                          // __imc_Coord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Coord_sys.hh
//---------------------------------------------------------------------------//
