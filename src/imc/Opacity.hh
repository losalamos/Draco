//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Opacity.hh
 * \author Thomas M. Evans
 * \date   Fri Feb  6 13:52:29 1998
 * \brief  Opacity class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Opacity_hh__
#define __imc_Opacity_hh__

#include "Frequency.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <iomanip>
#include <iostream>

namespace rtt_imc 
{

//===========================================================================//
/*!

 * \class Opacity

 * \brief Opacity class for Fleck and Cumming's IMC.

 * The Opacity class stores the opacities necessary for one timestep of
 * transport.  It is designed to be rebuilt each timestep when the material
 * properties cause a change in the particle opacities, although this is not
 * required.  The data it provides are:

 * \arg absorption opacities (1/cm)
 * \arg coherent (isotropic) scattering cross sections (1/cm)
 * \arg fleck factor (dimensionless)

 * The Opacity class simply stores these opacities for IMC.  It is the job of
 * a builder class to actually generate the opacities that Opacity stores.
 * Thus, for gray calculations, the Planck and absorption opacities have the
 * same value.

 * Additionally, the Opacity class uses the stored opacities to return Fleck
 * and Cummings specific "opacity-like" quantities. These result from the
 * time-implicit differencing in the material-energy equation and the
 * time-continuous treatment of the transport equation.  These are:

 * \arg the effective absorption
 * \arg the effective scattering
 */
/*!
 * \example imc/test/tstOpacity.cc

 * Example usage of the rtt_imc::Opacity class.

 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class FT>
class Opacity
{
  private:
    Opacity();
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

template<class MT, class FT>
std::ostream& operator<<(std::ostream &output, const Opacity<MT,FT> &object)
{
    object.print(output);
    return output;
}

//===========================================================================//
/*!
 * \brief Specialization of Opacity class for Gray_Frequency.
 */
//===========================================================================//

template<class MT>
class Opacity<MT, Gray_Frequency>
{
  public:
    // Useful typedefs.
    typedef typename MT::template CCSF<double>    ccsf_double;
    typedef rtt_dsxx::SP<Gray_Frequency>          SP_Frequency;

  private:
    // Absorption opacities in 1/cm.
    ccsf_double sigma_abs;

    // Scattering (coherent) cross sections in 1/cm.
    ccsf_double sigma_thomson;

    // Fleck factor.
    ccsf_double fleck;

    // Frequency definitions.
    SP_Frequency frequency;
    
  public:
    // Opacity constructor.
    Opacity(SP_Frequency, const ccsf_double &, const ccsf_double &,
	    const ccsf_double &);

    // >>> ACCESSORS

    // Get the group absorption opacity in a cell (/cm).
    inline double get_sigma_abs(int cell) const;

    // Get the group (coherent) scattering cross section in a cell (/cm).  
    inline double get_sigma_thomson(int cell) const;

    //! Get the Fleck factor in a cell.
    double get_fleck(int cell) const { return fleck(cell); }

    //! Get the number of cells that these opacities are stored in.
    int num_cells() const { return sigma_abs.get_Mesh().num_cells(); }

    //! Get a Smart Pointer to the frequency.
    SP_Frequency get_Frequency() const { return frequency; }

    // Get integrated normalized Planck function in a cell (unitless).
    double get_integrated_norm_Planck(int cell) const { return 1.0; }

    // >>> FLECK AND CUMMINGS OPACITY OPERATIONS

    // Get the effective scattering cross section.
    inline double get_sigeffscat(int cell) const;

    // Get the effective absorption cross section.
    inline double get_sigeffabs(int cell) const;

    // Print diagnostic function.
    void print(std::ostream &) const;
};

//---------------------------------------------------------------------------//
// FLECK AND CUMMINGS OPACITY OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return a cell-centered gray absorption opacity in per cm.
 *
 * \param cell cell index [1,N_cells]
 */
template<class MT>
double Opacity<MT,Gray_Frequency>::get_sigma_abs(int cell) const
{
    Require (cell > 0 && cell <= num_cells());
    return sigma_abs(cell);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a cell-centered gray scattering cross section in per cm.
 *
 * \param cell cell index [1,N_cells]
 */
template<class MT>
double Opacity<MT,Gray_Frequency>::get_sigma_thomson(int cell) 
    const
{
    Require (cell > 0 && cell <= num_cells());
    return sigma_thomson(cell);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the Fleck and Cummings effective scattering cross section.

 * The effective scattering cross section in the Fleck and Cummings IMC
 * method is sigma = (1-F)*sigma_absorption, where F is the Fleck factor.
 *
 * \param cell cell index [1,N_cells]

 */
template<class MT>
double Opacity<MT,Gray_Frequency>::get_sigeffscat(int cell) const 
{
    Require (cell > 0 && cell <= num_cells());
    return (1.0 - fleck(cell)) * sigma_abs(cell);
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Return the Fleck and Cummings effective absorption cross section.
  
 * The effective absorption cross section in the Fleck and Cummings IMC
 * method is sigma = F*sigma_absorption, where F is the Fleck factor.
 *
 * \param cell cell index [1,N_cells]
 
 */
template<class MT>
double Opacity<MT,Gray_Frequency>::get_sigeffabs(int cell) const 
{ 
    Require (cell > 0 && cell <= num_cells());
    return fleck(cell) * sigma_abs(cell);
}

//===========================================================================//
/*!
 * \brief Specialization of Opacity class for Multigroup_Frequency.
 */
//===========================================================================//

template<class MT>
class Opacity<MT, Multigroup_Frequency>
{
  public:
    // Useful typedefs.
    typedef std::vector<double>                   sf_double;
    typedef std::vector<sf_double>                vf_double;
    typedef typename MT::template CCSF<sf_double> ccsf_vector;
    typedef typename MT::template CCSF<double>    ccsf_double;
    typedef rtt_dsxx::SP<Multigroup_Frequency>    SP_Frequency;

  private:
    // Absorption opacities in 1/cm.
    ccsf_vector sigma_abs;

    // Scattering (coherent) cross sections in 1/cm.
    ccsf_vector sigma_thomson;

    // Fleck factor.
    ccsf_double fleck;

    // Integrated normalized Planck function per cell.
    ccsf_double integrated_norm_planck;

    // Unnormalized cdf of sigma * normalized Planck
    ccsf_vector emission_group_cdf;

    // Frequency definitions.
    SP_Frequency frequency;
    
  public:
    // Opacity constructor.
    Opacity(SP_Frequency, const ccsf_vector &, const ccsf_vector &, 
	    const ccsf_double &, const ccsf_double &, const ccsf_vector &);

    // >>> ACCESSORS

    // Get the group absorption opacity in a cell (/cm).
    inline double get_sigma_abs(int cell, int grp) const;

    // Get the group (coherent) scattering cross section in a cell (/cm).  
    inline double get_sigma_thomson(int cell, int grp) const;

    //! Get the Fleck factor in a cell.
    double get_fleck(int cell) const { return fleck(cell); }

    //! Get the number of cells that these opacities are stored in.
    int num_cells() const { return sigma_abs.get_Mesh().num_cells(); }

    //! Get a Smart Pointer to the frequency.
    SP_Frequency get_Frequency() const { return frequency; }

    // Get integrated normalized Planck function in a cell (unitless).
    inline double get_integrated_norm_Planck(int cell) const;

    // Get emission group Cumulative Distribution Function.
    inline sf_double get_emission_group_cdf(int c) const;

    // >>> FLECK AND CUMMINGS OPACITY OPERATIONS

    // Get the effective scattering cross section.
    inline double get_sigeffscat(int cell, int grp) const;

    // Get the effective absorption cross section.
    inline double get_sigeffabs(int cell, int grp) const;

    // Print diagnostic function.
    void print(std::ostream &) const;
};

//---------------------------------------------------------------------------//
// FLECK AND CUMMINGS OPACITY OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return a cell-centered group absorption opacity in per cm.
 *
 * \param cell cell index [1,N_cells]
 * \param group group index [1,N_groups]
 */
template<class MT>
double Opacity<MT, Multigroup_Frequency>::get_sigma_abs(int cell, int group)
    const
{
    Require (cell > 0 && cell <= num_cells());
    
    Require (group > 0 && group <= frequency->get_num_groups());
    Check (frequency->get_num_groups() == sigma_abs(cell).size());
    return sigma_abs(cell)[group-1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return a cell-centered group scattering cross section in per cm.
 *
 * \param cell cell index [1,N_cells]
 * \param group group index [1,N_groups]
 */
template<class MT>
double Opacity<MT,Multigroup_Frequency>::get_sigma_thomson(int cell, 
							   int group) const
{
    Require (cell > 0 && cell <= num_cells());
    Require (group > 0 && group <= frequency->get_num_groups());
    Check (frequency->get_num_groups() == sigma_abs(cell).size());

    return sigma_thomson(cell)[group-1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the Fleck and Cummings effective scattering cross section.

 * The effective scattering cross section in the Fleck and Cummings IMC
 * method is sigma = (1-F)*sigma_absorption, where F is the Fleck factor.
 *
 * \param cell cell index [1,N_cells]
 * \param group group index [1,N_groups]

 */
template<class MT>
double Opacity<MT,Multigroup_Frequency>::get_sigeffscat(int cell,
							int group) const 
{
    Require (cell > 0 && cell <= num_cells());
    Require (group > 0 && group <= frequency->get_num_groups());
    Check (frequency->get_num_groups() == sigma_abs(cell).size());

    return (1.0 - fleck(cell)) * sigma_abs(cell)[group-1];
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Return the Fleck and Cummings effective absorption cross section.
  
 * The effective absorption cross section in the Fleck and Cummings IMC
 * method is sigma = F*sigma_absorption, where F is the Fleck factor.
 *
 * \param cell cell index [1,N_cells]
 * \param group group index [1,N_groups]
 
 */
template<class MT>
double Opacity<MT,Multigroup_Frequency>::get_sigeffabs(int cell, int group)
    const 
{ 
    Require (cell > 0 && cell <= num_cells());
    Require (group > 0 && group <= frequency->get_num_groups());
    Check (frequency->get_num_groups() == sigma_abs(cell).size());

    return fleck(cell) * sigma_abs(cell)[group-1];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get integrated normalized Planck function in a cell (unitless).
 */
template<class MT>
double Opacity<MT,Multigroup_Frequency>::get_integrated_norm_Planck(
    int cell) const
{
    Require (cell > 0 && cell <= num_cells());

    return integrated_norm_planck(cell);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get emission group Cumulative Distribution Function.
 */
template<class MT>
typename Opacity<MT,Multigroup_Frequency>::sf_double 
Opacity<MT,Multigroup_Frequency>::get_emission_group_cdf(int cell) const
{
    Require (cell > 0 && cell <= num_cells());

    return emission_group_cdf(cell);
}

} // end namespace rtt_imc

#endif                          // __imc_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity.hh
//---------------------------------------------------------------------------//
