//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Coord_sys.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 30 16:45:37 1998
 * \brief  Coord_sys base class implementation file
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Coord_sys.hh"
#include "Math.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

namespace rtt_mc 
{

//---------------------------------------------------------------------------//
// VIRTUAL MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a direction on the unit sphere.
 *
 * This function samples a direction on the unit sphere for the following
 * distributions:
 * - "iostropic"
 * .
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
 * \param dist std::string for the directional distribution; presently only
 * "isotropic" is supported
 * \param random rtt_rng::Sprng random number object
 * \return vector<double> of [omega_x, omega_y, omega_z] direction cosines 
 */
Coord_sys::sf_double Coord_sys::sample_dir(std_string       dist, 
					   const rng_Sprng &random) const
{

    using rtt_dsxx::soft_equiv;
    using std::vector;
    using std::sqrt;

    // make return vector
    vector<double> omega_;

    // get direction cosines for different distributions
    if (dist == "isotropic")
    {
	omega_ = sample_isotropic_dir(random);
	Check (omega_.size() == 3);
    }
    else
    {
	Insist(0, "Only an isotropic distribution is available!");
    }

    // calculate the magnitude of the unit direction vector
    double magnitude = sqrt(rtt_mc::global::dot(omega_, omega_));

    // weak check on the magnitude
    Check (soft_equiv(magnitude, 1.0, 1.0e-5));

    // renormalize if magnitude fails stronger check
    if (!soft_equiv(magnitude, 1.0, 1.0e-10))
    {
	omega_[0] /= magnitude;
	omega_[1] /= magnitude;
	omega_[2] /= magnitude;
	
	double new_magnitude = sqrt(rtt_mc::global::dot(omega_, omega_));
	Check (soft_equiv(new_magnitude, 1.0, 1.0e-10));
    }

    // return vector
    return omega_;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Calculate new direction cosines when scattered through an angle
 * (theta, phi).
 *
 * The direction cosines are given by the formula:
 * \f[
 * \mathbf{\Omega} = \Omega_{x}\mathbf{i} + \Omega_{y}\mathbf{j} 
 * + \Omega_{z}\mathbf{k}
 * \f]
 * When scattered through an angle \f$(\theta,\phi)\f$, the new direction
 * cosines are:
 * \f[
 * \Omega_{x}' = \Omega_{x}\cos\theta + \Omega_{x}\Omega_{z}\sin\theta
 * \cos\phi / \alpha - \Omega_{y}\sin\theta\sin\phi/\alpha
 * \f]
 * \f[
 * \Omega_{y}' = \Omega_{z}\cos\theta + \Omega_{y}\Omega_{z}\sin\theta
 * \cos\phi / \alpha + \Omega_{x}\sin\theta\sin\phi/\alpha
 * \f]
 * \f[
 * \Omega_{z}' = \Omega_{z}\cos\theta - \alpha\sin\theta\cos\phi
 * \f]
 * where
 * \f[
 * \alpha = \sqrt{1-\Omega_{z}^{2}}
 * \f]
 * 
 * \param costheta cosine of polar scattering angle
 * \param phi azimuthal scattering angle
 * \param omega_ direction cosines on input, new direction cosines on output
 */
void Coord_sys::calc_omega(double     costheta, 
			   double     phi, 
			   sf_double &omega_) const
{
    using rtt_dsxx::soft_equiv;
    using std::vector;
    using std::sqrt;
    using std::cos;
    using std::sin;

    Require (omega_.size() == 3);

    // check that the old, incoming direction is reasonably normalized
    double old_magnitude = sqrt(rtt_mc::global::dot(omega_, omega_));
    Check (soft_equiv(old_magnitude, 1.0, 1.e-5));

    // renormalize if magnitude fails stronger check
    if (!soft_equiv(old_magnitude, 1.0, 1.0e-10))
    {
	omega_[0] /= old_magnitude;
	omega_[1] /= old_magnitude;
	omega_[2] /= old_magnitude;
	
	double new_magnitude = sqrt(rtt_mc::global::dot(omega_, omega_));
	Check (soft_equiv(new_magnitude, 1.0, 1.0e-10));
    }

    // calculate new direction cosines
    double sintheta        = sqrt(1.0 - costheta * costheta);
    vector<double> old_dir = omega_;
    double factor          = sqrt(std::fabs(1.0 - old_dir[2] * old_dir[2]));

    // if factor is zero use degenerate forms of scattering transformation
    if (factor < 1.e-6)
    {
        omega_[0] = sintheta * cos(phi);
        omega_[1] = sintheta * sin(phi);
        omega_[2] = old_dir[2] / std::fabs(old_dir[2]) * costheta;
    }
    
    // otherwise do standard transformation
    else
    {
        omega_[0] = old_dir[0] * costheta + old_dir[2] * old_dir[0] *
            sintheta *cos(phi) / factor - old_dir[1] * sintheta *
            sin(phi) / factor;
        omega_[1] = old_dir[1] * costheta + old_dir[2] * old_dir[1] *
            sintheta *cos(phi) / factor + old_dir[0] * sintheta *
            sin(phi) / factor;
        omega_[2] = old_dir[2] * costheta - factor * sintheta * cos(phi);
    }

    // calculate the magnitude of the new unit direction vector
    double magnitude = sqrt(rtt_mc::global::dot(omega_, omega_));

    // weak check on the magnitude
    Check (soft_equiv(magnitude, 1.0, 1.0e-5));

    // renormalize if magnitude fails stronger check
    if (!soft_equiv(magnitude, 1.0, 1.0e-10))
    {
	omega_[0] /= magnitude;
	omega_[1] /= magnitude;
	omega_[2] /= magnitude;
	
	double new_magnitude = sqrt(global::dot(omega_, omega_));
	Check (soft_equiv(new_magnitude, 1.0, 1.0e-10));
    }
}

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Coord_sys.cc
//---------------------------------------------------------------------------//
