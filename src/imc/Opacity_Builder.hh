//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.hh
// Thomas M. Evans
// Fri Mar  6 17:21:36 1998
//---------------------------------------------------------------------------//
// @> Opacity_Builder class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Opacity_Builder_hh__
#define __imc_Opacity_Builder_hh__

//===========================================================================//
// class Opacity_Builder - 
//
// Purpose : build the Opacity and Mat_State objects, templated on MT
//
// revision history:
// -----------------
//  0) original
//  1)  7-28-98 : converted density, kappa, temperature, and specific_heat to 
//                cell-based arrays in the Interface, ie. we no longer need
//                zone or mat_zone
//  2)  7-31-01 : added analytic opacity offset, kappa_offset (tju)
// 
//===========================================================================//

#include "Opacity.hh"
#include "Mat_State.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <string>

namespace rtt_imc 
{

template<class MT>
class Opacity_Builder
{
    // typedefs
    typedef rtt_dsxx::SP<MT>             SP_MT;
    typedef rtt_dsxx::SP<Mat_State<MT> > SP_Mat;
    typedef rtt_dsxx::SP<Opacity<MT> >   SP_Opacity; 

  private:
    // data received from XX_Interface
    std::vector<double> density;
    std::vector<double> kappa;
    std::vector<double> kappa_offset;
    std::vector<double> kappa_thomson;
    std::vector<double> temperature;
    std::vector<double> specific_heat;
    std::vector<int>    material_id;
    double implicitness;
    double delta_t;
    std::string analytic_opacity;
    std::string analytic_sp_heat;

  public:
    // templated explicit constructor depends on interface type (IT)
    template<class IT>
    explicit Opacity_Builder(rtt_dsxx::SP<IT>);

    // build state member functions

    // build Mat_State helper functions
    SP_Mat build_Mat(SP_MT);
    
    // build Opacity helper functions
    SP_Opacity build_Opacity(SP_MT, SP_Mat);
};
    
//---------------------------------------------------------------------------//
// TEMPLATE MEMBERS
//---------------------------------------------------------------------------//
//! Contructor.

template<class MT>
template<class IT>
Opacity_Builder<MT>::Opacity_Builder(rtt_dsxx::SP<IT> interface)
{
    Require (interface);

    // assign data members from the interface parser
    density          = interface->get_density();
    kappa            = interface->get_kappa();
    kappa_offset     = interface->get_kappa_offset();
    kappa_thomson    = interface->get_kappa_thomson();
    temperature      = interface->get_temperature();
    specific_heat    = interface->get_specific_heat();
    implicitness     = interface->get_implicit();	
    delta_t          = interface->get_delta_t();
    analytic_opacity = interface->get_analytic_opacity();
    analytic_sp_heat = interface->get_analytic_sp_heat();

    // some crucial checks about our data
    Check (implicitness >= 0 && implicitness <= 1);
    Check (delta_t > 0);
}

} // end namespace rtt_imc

#endif                          // __imc_Opacity_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity_Builder.hh
//---------------------------------------------------------------------------//
