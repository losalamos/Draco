//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  Pseudo_Line_Analytic_Odfmg_Opacity class definition.
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Pseudo_Line_Analytic_Odfmg_Opacity_hh__
#define __cdi_analytic_Pseudo_Line_Analytic_Odfmg_Opacity_hh__

#include "Analytic_Models.hh"
#include "Analytic_Odfmg_Opacity.hh"
#include "Pseudo_Line_Base.hh"
#include "cdi/OpacityCommon.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_cdi_analytic
{

//===========================================================================//
/*!
 * \class Pseudo_Line_Analytic_Odfmg_Opacity
 *
 * \brief Derived rtt_cdi::OdfmgOpacity class for analytic opacities.
 *
 * Primarily code from Analytic_Multigroup_Opacity.
 */
// 
//===========================================================================//

class Pseudo_Line_Analytic_Odfmg_Opacity :
        public Analytic_Odfmg_Opacity, public Pseudo_Line_Base
{
  private:

    Averaging averaging_;
    unsigned qpoints_;

  public:
    // Constructor.
    Pseudo_Line_Analytic_Odfmg_Opacity(
        const sf_double         &groups,
        const sf_double         &bands,
        rtt_cdi::Reaction        reaction_in,
        SP<Expression const> const &cont,
        unsigned number_of_lines,
        double line_peak,
        double line_width,
        unsigned number_of_edges,
        double edge_ratio,
        double Tref,
        double Tpow,
        double emin,
        double emax,
        Averaging averaging,
        unsigned qpoints,
        unsigned seed);

    // Constructor for packed Pseudo_Line_Analytic_Odfmg_Opacities
    explicit Pseudo_Line_Analytic_Odfmg_Opacity(const sf_char &);

    std::vector< std::vector<double> > getOpacity( 
        double targetTemperature,
        double targetDensity ) const; 

    std::vector< std::vector< std::vector<double> > > getOpacity( 
        const std::vector<double>& targetTemperature,
        double targetDensity ) const; 

    std::vector< std::vector< std::vector<double> > > getOpacity( 
        double targetTemperature,
        const std::vector<double>& targetDensity ) const; 

    // Get the data description of the opacity.
    inline std_string getDataDescriptor() const;

    // Pack the Pseudo_Line_Analytic_Odfmg_Opacity into a character string.
    sf_char pack() const;
};

} // end namespace rtt_cdi_analytic

#endif              // __cdi_analytic_Pseudo_Line_Analytic_Odfmg_Opacity_hh__

//---------------------------------------------------------------------------//
//                   end of cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.hh
//---------------------------------------------------------------------------//
