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

#include "Constants.hh"
#include "rng/Sprng.hh"
#include <vector>
#include <string>
#include <cmath>

namespace rtt_mc 
{

//===========================================================================/
/*!
 * \class Coord_sys
 *
 * Coordinate system abstract base class for use in mesh types for Monte
 * Carlo applications.  The coordinate system family of classes provide a
 * means for doing coordinate system dependent sampling and particle
 * direction calculations.  In general, the coordinate system classes are
 * used as components of mesh types; thus, they are not generally used
 * outside of mesh classes.  For examples of this type of implementation see
 * the rtt_mc::OS_Mesh class.
 *
 * One place the coordinate system is used outside of the MT (Mesh Type)
 * implementations is in rtt_imc::Random_Walk.  The Coord_sys is used to
 * sample a position on the surface of a sphere.
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
//  6) 5-SEP-01 : made checks consistent on magnitude of direction cosines.
//                For soft_equiv(*,1.0,eps), we use eps=1.0e-5 for
//                Omega-dot-Omega = magnitude^2 and eps=1.0e-10 for
//                sqrt(Omega-dot-Omega) = magnitude, since magnitude^2 is
//                farther away from unity than  magnitude.  Checks fail if
//                magnitude is more than 1.0e-5 from  unity.  Omega is
//                renormalized if its magnitude is more than 1.0e-10 from
//                unity.
//  7) 02-14-03 : added sampling on the surface of a sphere
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
    const int eff_dimension;

  public:
    //! Constructor.
    Coord_sys(int dimension_) : dimension(dimension_), eff_dimension(3) {}

    //! Virtual destructor to insure correct behavior down inheritance chain.
    virtual ~Coord_sys() {}

    // >>> BASE CLASS MEMBERS

    //! Return the dimensionality of the coordinate system.
    int get_dim() const { return dimension; } 

    //! Return the effective dimension of the system.
    int get_sdim() const { return eff_dimension; }

    // >>> COORDINATE SYSTEM INTERFACE 

    //! Return a string descriptor of the coordinate system.
    virtual std_string get_Coord() const = 0;
 
    //! Sample position uniformly in bounded 3-D space.
    virtual sf_double sample_pos(const sf_double &, const sf_double &, 
				 const rng_Sprng &) const = 0; 

    //! Sample position in bounded 3-D space in a linear distribution. 
    virtual sf_double sample_pos(const sf_double &, const sf_double &, 
				 const rng_Sprng &, const sf_double &, 
				 double) const = 0;

    //! Sample a position uniformly on a face.
    virtual sf_double sample_pos_on_face(const sf_double &, const sf_double &,
					 int, const rng_Sprng &) const = 0;

    // Sample a particle direction.
    virtual sf_double sample_dir(std_string, const rng_Sprng &) const;

    // Calculate a new particle direction for scattering.
    virtual void calc_omega(double, double, sf_double &) const;

    // Sample an isotropic particle direction.
    inline virtual sf_double sample_isotropic_dir(const rng_Sprng &) const;

    // >>> OVERLOADED OPERATORS

    // Equality operator.
    inline bool operator==(const Coord_sys &) const;

    //! Inequality operator.
    bool operator!=(const Coord_sys &rhs) const { return !(*this == rhs); }
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Sample an isotropic direction on the unit sphere.
 *
 * This function samples an isotropic direction on the unit sphere. It is
 * used by sample_dir for isotropic distributions.  sample_dir() contains
 * more checking on the magnitude of the direction cosines.
 *
 * The direction cosines are given by the formula:
 * \f[
 * \mathbf{\Omega} = \Omega_{x}\mathbf{i} + \Omega_{y}\mathbf{j} 
 * + \Omega_{z}\mathbf{k}
 * \f]
 * where
 * \f[
 * \Omega_{x} = \cos\theta\sin\phi
 * \f]
 * \f[
 * \Omega_{y} = \sin\theta\sin\phi
 * \f]
 * \f[
 * \Omega_{z} = \cos\theta
 * \f]
 * We sample \f$(\theta,\phi)\f$ according to the specified distribution and
 * return a std::vector<double> of \f$(\Omega_{x},\Omega_{y},\Omega_{z})\f$. 
 *
 * \param random rtt_rng::Sprng random number object
 * \return vector<double> of [omega_x, omega_y, omega_z] direction cosines 
 */
Coord_sys::sf_double Coord_sys::sample_isotropic_dir(const rng_Sprng &random) 
    const
{
    using std::vector;
    using std::sqrt;
    using std::sin;
    using std::cos;

    // make return vector
    vector<double> omega_(3);

    // sample costheta and phi for 3D transport 
    double costheta = 1.0 - 2.0 * random.ran();
    double sintheta = sqrt(1.0 - costheta * costheta);
    double phi      = 2.0 * global::pi * random.ran();

    // calculate 3D direction cosines
    omega_[0] = sintheta * cos(phi);
    omega_[1] = sintheta * sin(phi);
    omega_[2] = costheta;

    // return the direction
    return omega_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Equality operator.
 */
bool Coord_sys::operator==(const Coord_sys &rhs) const
{
    if (dimension == rhs.dimension && eff_dimension == rhs.eff_dimension)
	return true;
    return false;
}

} // end namespace rtt_mc

#endif                          // __mc_Coord_sys_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Coord_sys.hh
//---------------------------------------------------------------------------//
