//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Diffusion_Opacity.hh
 * \author Thomas M. Evans
 * \date   Wed Feb 19 10:43:46 2003
 * \brief  Diffusion_Opacity class definition
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef RTT_imc_Diffusion_Opacity_t_HH
#define RTT_imc_Diffusion_Opacity_t_HH

#include "Fleck_Factors.hh"
#include "mc/Constants.hh"
#include "ds++/SP.hh"

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Diffusion_Opacity
 *
 * \brief Diffusion opacities used in hybrid diffusion/IMC methods.
 *
 * The Diffusion_Opacity class holds diffusion opacity data that is used by
 * hybrid diffusion/IMC methods.  Chief among these methods are random walk
 * (rtt_imc::Random_Walk) and discrete diffusion (DDIMC).  The data held by
 * the class are Fleck factors and gray Rosseland opacities (needed to
 * generate diffusion coefficients):
 * 
 * - Rosseland gray absorption opacities (1/cm)
 * - fleck factors (dimensionaless)
 *
 * Also, the class can generate diffusion coefficients and the gray Rosseland
 * effective scattering opacity (1/cm), (1-fleck)$\sigma_R$, for Random Walk
 * and, in the future, DDIMC.
 *
 * The Diffusion_Opacity class is not frequency-dependent.  All diffusion
 * data in this class is "gray".
 *
 * The Diffusion_Opacity class is templated on mesh-type (MT).
 *
 */
// revision history:
// -----------------
// 0) original
// 1) 25-AUG-03 : updated to hold the Rossealand absorption opacity (instead
//                of Rosseland total opacity) and the total gray scattering
//                opacity
// 
//===========================================================================//

template<class MT>
class Diffusion_Opacity 
{
  public:
    // Useful typedefs.
    typedef typename MT::template CCSF<double> ccsf_double;
    typedef rtt_dsxx::SP<Fleck_Factors<MT> >   SP_Fleck_Factors;

  private:
    // >>> DATA
    
    // Fleck factors.
    SP_Fleck_Factors fleck;

    // Rosseland gray absorption opacity.
    ccsf_double rosseland;

    // Total gray scattering opacity.
    ccsf_double scattering;

  public:
    // Constructor.
    Diffusion_Opacity(SP_Fleck_Factors, const ccsf_double &, 
		      const ccsf_double &);

    // >>> ACCESSORS

    //! Get the number of cells that these opacities are stored in.
    int num_cells() const { return rosseland.get_Mesh().num_cells(); }
    
    //! Get the Fleck factor in a cell.
    double get_fleck(int cell) const { return fleck->fleck(cell); }

    //! Get the gray Rosseland absorption opacity in a cell (1/cm).
    double get_Rosseland_opacity(int c) const { return rosseland(c); }

    //! Get the total, gray scattering opacity in a cell (1/cm).
    double get_gray_scattering(int c) const { return scattering(c); }

    // Get the effective Rosseland mean free path (cm).
    inline double get_Rosseland_effmfp(int cell) const;

    // Get diffusion coefficents for random walk per cell.
    inline double get_random_walk_D(int cell) const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the random walk diffusion coefficent per cell in
 * (cm^2/shake). 
 *
 * The random walk approximation utilizes the following diffusion coefficient
 * (Fleck and Canfield, 1984):
 * \f[
 * D = \frac{c}{3(1-f)\sigma_{R}}
 * \f]
 * where \f$ c \f$ is the speed of light in [cm/shake], \f$ f \f$ is the
 * Fleck factor, and \f$\sigma_{R}\f$ is the Rosseland gray opacity in
 * [1/cm].  Thus, \f$ D \f$ has units of [cm^2/shake].
 */
template<class MT>
double Diffusion_Opacity<MT>::get_random_walk_D(int cell) const
{
    using rtt_mc::global::c;

    Require (cell > 0 && cell <= num_cells());
    Check (fleck->fleck(cell) >= 0.0);
    Check (fleck->fleck(cell) <= 1.0);

    return c / (3.0 * ((1.0-fleck->fleck(cell)) * rosseland(cell) + 
		scattering(cell))); 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the effective Rosseland mean free path.
 *
 * The effective Rosseland mean free path is:
 * \f[
 * \lambda_{R\mbox{eff}} = \frac{1}{(1-f)\sigma_{R\mbox{abs}} + 
 * \sigma_{s}}
 * \f]
 * where the scattering cross section includes both Thomson and Compton
 * scattering. 
 */
template<class MT>
double Diffusion_Opacity<MT>::get_Rosseland_effmfp(int cell) const
{

    Require (cell > 0 && cell <= num_cells());
    Check (fleck->fleck(cell) >= 0.0);
    Check (fleck->fleck(cell) <= 1.0);

    double denominator = (1.0-fleck->fleck(cell)) * rosseland(cell) +
	scattering(cell);

    // make sure that we don't divide by zero
    if (denominator == 0.0)
	return rtt_mc::global::huge;

    return 1.0 / denominator;
}

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * \param fleck_in rtt_dsxx::SP to a rtt_imc::Fleck_Factors object
 * \param rosse_in MT::CCSF of gray Rosseland opacities in (1/cm)
 */
template<class MT>
Diffusion_Opacity<MT>::Diffusion_Opacity(SP_Fleck_Factors   fleck_in,
					 const ccsf_double &rosse_in,
					 const ccsf_double &scat_in)
    : fleck(fleck_in),
      rosseland(rosse_in),
      scattering(scat_in)
{
    Require (fleck);
    Require (fleck->fleck.size() == rosseland.get_Mesh().num_cells());
    Require (fleck->fleck.size() == scattering.get_Mesh().num_cells());
}

} // end namespace rtt_imc

#endif                // RTT_imc_Diffusion_Opacity_t_HH

//---------------------------------------------------------------------------//
//                              end of imc/Diffusion_Opacity.hh
//---------------------------------------------------------------------------//
