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
// 
//===========================================================================//

#include "Opacity.hh"
#include "Mat_State.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_imc 
{

template<class MT>
class Opacity_Builder
{
    // typedefs
    typedef dsxx::SP<MT> SP_MT;
    typedef dsxx::SP<Mat_State<MT> > SP_Mat;
    typedef dsxx::SP<Opacity<MT> > SP_Opacity; 

  private:

    // data received from XX_Interface
    std::vector<double> density;
    std::vector<double> kappa;
    std::vector<double> kappa_thomson;
    std::vector<double> temperature;
    std::vector<double> specific_heat;
    std::vector<int>    material_id;
    double implicitness;
    double delta_t;
    std::string analytic_opacity;
    std::string analytic_sp_heat;

    // Begin_Doc opacity_builder-int.tex
    // Begin_Verbatim 

  public:
    // templated explicit constructor depends on interface type (IT)
    template<class IT>
    explicit Opacity_Builder(dsxx::SP<IT>);

    // build state member functions

    // build Mat_State helper functions
    SP_Mat build_Mat(SP_MT);
    
    // build Opacity helper functions
    SP_Opacity build_Opacity(SP_MT, SP_Mat);	

    // End_Verbatim 
    // End_Doc 
};
    
} // end namespace rtt_imc

#endif                          // __imc_Opacity_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Opacity_Builder.hh
//---------------------------------------------------------------------------//
