//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/XYZCoord_sys.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 30 16:45:36 1998
 * \brief  XYZCoord_sys derived class header file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_XYZCoord_sys_hh__
#define __mc_XYZCoord_sys_hh__

#include "Coord_sys.hh"
#include "Constants.hh"
#include "rng/Sprng.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <vector>
#include <string>
#include <cmath>

namespace rtt_mc 
{

//===========================================================================//
/*!
 * \class XYZCoord_sys
 *
 * Derived class for \i XYZ coordinate systems.  Inherits the public
 * interface from rtt_mc::Coord_sys.  
 */
// revision history:
// -----------------
//  0) original
//  1)  2-18-98 : added getCoord() function for debugging
//  2)  3-12-98 : moved Calc and Set_omega functions to Coord_sys because
//                they are the same for XY and XYZ transport
//  3)  3-16-98 : reserve calc_normal function for later if need be
//  4)   5-5-98 : added sample_pos virtual function
//  5)  6-10-98 : added sample_pos_on_face virtual function
//  6)  6-12-98 : changed interface to sample_pos()
//  7)  4-13-99 : moved to mc package
//  8)  2-14-03 : added sample_pos_on_sphere virtual function to class
// 
//===========================================================================//

class XYZCoord_sys : public Coord_sys
{
  public:
    // Default constructor for 3D meshes.
    XYZCoord_sys() : Coord_sys(3) {}

    // >>> VIRTUAL FUNCTIONS INHERITED FROM COORD_SYS

    //! Return the coordinate system.
    std_string get_Coord() const { std_string c = "xyz"; return c; }

    // Sample position uniformly in XYZ coordinate system.
    inline sf_double sample_pos(const sf_double &, const sf_double &, 
				const rng_Sprng &) const; 
    
    // Sample position in a linear distribution in XYZ coordinate system.
    inline sf_double sample_pos(const sf_double &, const sf_double &, 
				const rng_Sprng &, const sf_double &,
				double) const;
    
    // Sample position uniformly on a face.
    inline sf_double sample_pos_on_face(const sf_double &, const sf_double &, 
					int, const rng_Sprng &) const; 
};

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE TO XYZ COORDINATE SYSTEM (inherited from Coord_sys)
//---------------------------------------------------------------------------//
/*!
 * \brief Sample position unformly in bounded \e XYZ space.
 *
 * This function samples a position uniformly in \e XYZ space from [min,max]
 * along each dimension.  It is equivalent to sampling within a
 * parrallelpiped with bounds (min[0], max[0], min[1], max[1], min[2],
 * max[2]).
 * 
 * \param min minimum boundaries in \e (XYZ)
 * \param max maximum boundaries in \e (XYZ)
 * \param random rtt_rng::Sprng random number object
 * \return vector<double> of \e (XYZ) position
 */ 
Coord_sys::sf_double XYZCoord_sys::sample_pos(const sf_double &min, 
					      const sf_double &max,
					      const rng_Sprng &random) const
{
    // make return vector
    sf_double r(3);

    // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);

    for (int d = 0; d < 3; d++)
    {
	// do uniform sampling
	r[d] = (max[d] - min[d]) * random.ran() + min[d];
    }

    // return assigned array
    return r;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample position in bounded \e XYZ space in a linear distribution.
 * 
 * This function samples a position according to a linear distribution in \e
 * XYZ space from [min,max] along each dimension.  It is equivalent to
 * sampling within a parrallelpiped with bounds (min[0], max[0], min[1],
 * max[1], min[2], max[2]).
 * 
 * \param min minimum boundaries in \e (XYZ)
 * \param max maximum boundaries in \e (XYZ)
 * \param random rtt_rng::Sprng random number object
 * \param slope slope of linear distribution in \e (XYZ)
 * \param center_pt center points of linear distribution in \e (XYZ)
 * \return vector<double> of \e (XYZ) position
 */
Coord_sys::sf_double XYZCoord_sys::sample_pos(const sf_double &min, 
					      const sf_double &max,
					      const rng_Sprng &random, 
					      const sf_double &slope, 
					      double           center_pt) const
{
    // make return vector
    sf_double r(3);

    // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);

    for (int d = 0; d < 3; d++)
    {
	// sample the linear function using linear-linear decomposition of
	// y = mx + b (b is intercept on low side of cell)
	double b = center_pt - slope[d] * (max[d] - min[d]) * 0.5;

	// prob is the fractional area of the negative slope line
	double prob = 0.5 * b / center_pt;

	// sample the dimension
	if (random.ran() <= prob)
	    r[d] = max[d] - (max[d] - min[d]) * std::sqrt(random.ran());
	else
	    r[d] = min[d] + (max[d] - min[d]) * std::sqrt(random.ran());
    }

    // return assigned array
    return r;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Sample a position uniformly on a face in \e XYZ space.
 *
 * This function should be deprecated (move its implementation to OS_Mesh).
 * 
 * \param face face index (low X = 1, high X = 2, low Y = 3, high Y = 4, 
 * low Z = 5, high Z = 6)
 * \return vector<double> of \e (XYZ) position
 */
Coord_sys::sf_double XYZCoord_sys::sample_pos_on_face(
    const sf_double &min, 
    const sf_double &max, 
    int              face, 
    const rng_Sprng &random) const
{
    // make return vector
    sf_double r(3);

    // some assertions
    Check (min.size() == 3);
    Check (max.size() == 3);
    Check (face >= 1 && face <= 6);

    // distribute uniformly over face
    if (face == 1)
    {
	r[0] = min[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    else if (face == 2)
    {
	r[0] = max[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    else if (face == 3)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = min[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    else if (face == 4)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = max[1];
	r[2] = (max[2] - min[2]) * random.ran() + min[2];
    }
    else if (face == 5)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = min[2];
    }
    else if (face == 6)
    {
	r[0] = (max[0] - min[0]) * random.ran() + min[0];
	r[1] = (max[1] - min[1]) * random.ran() + min[1];
	r[2] = max[2];
    }

    // return assigned array
    return r;
}

} // end namespace rtt_mc

#endif                          // __mc_XYZCoord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of mc/XYZCoord_sys.hh
//---------------------------------------------------------------------------//
