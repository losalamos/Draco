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

#include "ds++/Assert.hh"
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
 * \arg planck (one-group) opacities (1/cm)
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
 * \arg the Fleck factor times the Planck absorption opacity

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

template<class MT>
class Opacity
{
  public:
    // Useful typedefs.
    typedef typename MT::CCSF_double ccsf_double;

  private:
    // Absorption opacity in 1/cm.
    ccsf_double sigma_abs;

    // Scattering (coherent) in 1/cm.
    ccsf_double sigma_thomson;

    // Plankian opacity in 1/cm.
    ccsf_double planck;

    // Fleck factor.
    ccsf_double fleck;
    
  public:
    // Opacity constructor.
    Opacity(const ccsf_double &, const ccsf_double &, const ccsf_double &, 
	    const ccsf_double &);

    // >>> ACCESSORS

    //! Get the absorption opacity in a cell in per cm.
    double get_sigma_abs(int cell) const { return sigma_abs(cell); }

    //! Get the (coherent) scattering cross section in a cell in per cm.
    double get_sigma_thomson(int cell) const { return sigma_thomson(cell); }

    //! Get the Planck absorption opacity in a cell in per cm.
    double get_planck(int cell) const { return planck(cell); }

    //! Get the Fleck factor in a cell.
    double get_fleck(int cell) const { return fleck(cell); }

    //! Get the number of cells that these opacities are stored in.
    int num_cells() const { return sigma_abs.get_Mesh().num_cells(); }

    // >>> FLECK AND CUMMINGS OPACITY OPERATIONS

    //! Return the Fleck factor times the Planck opacity in a cell.
    double fplanck(int cell) const { return fleck(cell) * planck(cell); }

    // Get the effective scattering cross section.
    inline double get_sigeffscat(int cell) const;

    // Get the effective absorption cross section.
    inline double get_sigeffabs(int cell) const;

    // Print diagnostic function.
    void print(std::ostream &, int) const;
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

template<class MT>
std::ostream& operator<<(std::ostream &output, const Opacity<MT> &object)
{
    // print out opacities for all cells
    using std::endl;
    using std::ios;
    using std::setw;
    using std::setiosflags;

    output << setw(8) << setiosflags(ios::right) << "Cell" 
	   << setw(15) << "Abs Opacity" 
	   << setw(15) << "Thomson Opac" << endl;
    output << "--------------------------------------" << endl;

    for (int i = 1; i <= object.num_cells(); i++)
        object.print(output, i);
    return output;
}

//---------------------------------------------------------------------------//
// FLECK AND CUMMINGS OPACITY OPERATIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the Fleck and Cummings effective scattering cross section.

 * The effective scattering cross section in the Fleck and Cummings IMC
 * method is sigma = (1-F)*sigma_absorption, where F is the Fleck factor.

 */
template<class MT>
double Opacity<MT>::get_sigeffscat(int cell) const 
{
    Require (cell > 0 && cell <= num_cells());
    return (1.0 - fleck(cell)) * sigma_abs(cell);
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Return the Fleck and Cummings effective absorption cross section.
  
 * The effective absorption cross section in the Fleck and Cummings IMC
 * method is sigma = F*sigma_absorption, where F is the Fleck factor.
 
 */
template<class MT>
double Opacity<MT>::get_sigeffabs(int cell) const 
{ 
    Require (cell > 0 && cell <= num_cells());
    return fleck(cell) * sigma_abs(cell);
}

} // end namespace rtt_imc

#endif                          // __imc_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity.hh
//---------------------------------------------------------------------------//
