//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Flat_Data_Container.hh
 * \author Thomas M. Evans
 * \date   Wed Feb  6 10:55:32 2002
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Flat_Data_Container_hh__
#define __imc_Flat_Data_Container_hh__

#include <vector>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Flat_Data_Container
 *
 * \brief Holds "flat" opacity and EOS data.
 *
 * This struct holds "flat" cell-centered fields of opacity and EOS data.  It
 * is used by the Flat_Mat_State_Builder and Flat_Data_Interface.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

struct Flat_Data_Container 
{
    // >>> TYPEDEFS
    typedef std::vector<double>    sf_double;
    typedef std::vector<sf_double> vf_double;
    
    // >>> DATA

    //! Cell-centered field of multigroup absorption opacities in /cm.
    vf_double mg_absorption_opacity;

    //! Cell-centered field of multigroup scattering opacities in /cm.
    vf_double mg_scattering_opacity;

    //! Group boundaries in keV.
    sf_double group_boundaries;

    //! Cell-centered field of gray absorption opacities in /cm.
    sf_double gray_absorption_opacity;

    //! Cell-centered field of gray scattering opacities in /cm.
    sf_double gray_scattering_opacity;

    //! Specific heats in Jerks/gm/keV.
    sf_double specific_heat;
};

} // end namespace rtt_imc

#endif                          // __imc_Flat_Data_Container_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Flat_Data_Container.hh
//---------------------------------------------------------------------------//
