//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Coord_sys.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 30 16:36:51 1998
 * \brief  Coord_sys abstract base class header file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Coord_sys_hh__
#define __mc_Coord_sys_hh__

#include "rng/Sprng.hh"
#include <vector>
#include <string>

namespace rtt_mc 
{

//===========================================================================/
/*!
 * \class Coord_sys
 *
 * Coordinate system abstract base class for use in mesh types for Monte
 * Carlo applications.  The coordinate system family of classes provide a
 * means for doing coordinate system dependent sampling.  In general, the
 * coordinate system classes are used as components of mesh types; thus, they
 * are not generally used outside of mesh classes.  For examples of this type
 * of implementation see the rtt_mc::OS_Mesh class.
 *
 * Currently, a majority of these functions will only work in orthogonal
 * meshes.  However, they can still be used for a whole class of orthogonal
 * meshes: \i XYZ, \i XY, \i RZ, \i R. 
 */
/*!
 * \example mc/test/tstCoord.cc
 *
 * Examples of coordinate system usage.  Normally these are used as
 * components of class types. 
 */
// revision history:
// -----------------
//  0) original
//  1)  3-12-98 : moved Calc and Set_omega functions into Coord_sys as
//                non-pure virtual functions because they are the same in
//                both XY and XYZ transport, added a transform for 2D meshes
//  2)  3-16-98 : reserve calc_normal function for later if need be
//  3)  3-17-98 : because of a dumb-ass oversight on my part, we don't need
//                a transform for 2D XY, it has been removed
//  4)  4-13-99 : moved into mc package
//  5) 3-AUG-00 : fixed calc_omega error! the degenerate z-case was incorrect
// 
//===========================================================================//

class Coord_sys
{
  public:
    // STL Typedefs
    typedef std::vector<int>    sf_int;
    typedef std::vector<double> sf_double;
    typedef std::string         std_string;
    typedef rtt_rng::Sprng      rng_Sprng;

  private:
    // Dimension of system.
    const int dimension;

    // Effective dimension of system (MC always tracks in 3-D).
    const int set_dimension;

  public:
    // Constructor for setting dimension of Coord_sys.
    Coord_sys(int dimension_) 
	:dimension(dimension_), set_dimension(3) {}

    // Virtual destructor to insure correct behavior down inheritance chain.
    virtual ~Coord_sys() {}

    // >>> Base class member functions.

    // We have two dimensionalities, a "real" dimension for the geometry and a
    // "transport" dimension for MC transport which is inherently 3D.
    int get_dim() const { return dimension; } 
    int get_sdim() const { return set_dimension; }

    // >>> Pure virtual functions.

    // Return the coordinate system.
    virtual std_string get_Coord() const = 0;
 
    // Sample positions based on the Coordinate system.
    virtual 
    sf_double sample_pos(sf_double &, sf_double &, rng_Sprng &) const = 0; 

    virtual 
    sf_double sample_pos(sf_double &, sf_double &, rng_Sprng &, sf_double &, 
			 double) const = 0;
    virtual 
    sf_double sample_pos_on_face(sf_double &, sf_double &, 
				 int, rng_Sprng &) const = 0;

    // >>> Virtual functions.

    // Sample directions.
    virtual sf_double sample_dir(std_string, rng_Sprng &) const;
    virtual void calc_omega(double, double, sf_double &) const;

    // Overloaded operators for equality.
    inline bool operator==(const Coord_sys &) const;
    bool operator!=(const Coord_sys &rhs) const { return !(*this == rhs); }
};

//---------------------------------------------------------------------------//
// OVERLOADED EQUALITY OPERATORS
//---------------------------------------------------------------------------//
// equality of Coordinate systems, because the derived classes contain no
// data we just have to worry about the base class type

bool Coord_sys::operator==(const Coord_sys &rhs) const
{
    if (dimension == rhs.dimension && set_dimension == rhs.set_dimension)
	return true;
    return false;
}

} // end namespace rtt_mc

#endif                          // __mc_Coord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Coord_sys.hh
//---------------------------------------------------------------------------//
