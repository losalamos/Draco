//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Flat_Mat_State_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:36:42 2001
 * \brief  Flat_Mat_State_Builder class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Flat_Mat_State_Builder_hh__
#define __imc_Flat_Mat_State_Builder_hh__

#include "Mat_State_Builder.hh"
#include <vector>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Flat_Mat_State_Builder
 *
 * This class builds instances of rtt_imc::Opacity and rtt_imc::Mat_State.
 * It receives its data from a defined interface in a "flat" form.  In other
 * words, the data fields are already defined.  It receives an Interface Type
 * (IT) in its constructor.  The IT must be a derived class of
 * rtt_imc::Interface and rtt_imc::Flat_Data_Interface.
 *
 * This class, along with rtt_imc::CDI_Mat_State_Builder, is a derived class
 * of Mat_State_Builder.  It should be used when a client already has
 * cell-centered material defined.  
 *
 * \sa rtt_imc::Interface, rtt_imc::Flat_Data_Interface.  See
 * tstMat_Data_Builder for examples of usage.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class Flat_Mat_State_Builder : public Mat_State_Builder<MT>
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>             SP_Mesh;
    typedef rtt_dsxx::SP<Mat_State<MT> > SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT> >   SP_Opacity;
    typedef std::vector<double>          sf_double;

  private:
    // Flat, cell-centered data fields received from the interface.

    // Densities in g/cc.
    sf_double density;

    // Absorption opacity in /cm.
    sf_double absorption_opacity;

    // Scattering opacity in /cm.
    sf_double scattering_opacity;

    // Material temperatures in keV.
    sf_double temperature;

    // Material specific heats in Jerks/gm/keV.
    sf_double specific_heat;

    // Fleck and Cummings implicitness factor.
    double    implicitness;

    // Timestep in shakes.
    double    delta_t;
    

  public:
    // Constructor.
    template<class IT>
    explicit Flat_Mat_State_Builder(rtt_dsxx::SP<IT>);

    // >>> PUBLIC INTERFACE

    // Build the Mat_State.
    SP_Mat_State build_Mat_State(SP_Mesh) const;

    // Build the Opacity.
    SP_Opacity build_Opacity(SP_Mesh m, SP_Mat_State s) const;
};

//---------------------------------------------------------------------------//
// MEMBER TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * The constructor gets data from the interface type.  The interface type
 * must provide member functions defined by rtt_imc::interface and
 * rtt_imc::Flat_Data_Interface.  However, because this is used as a template
 * argument, it is not required that the interface type inherit from these
 * classes.  It is only required that the interface type have these functions
 * defined.
 *
 * \param interface rtt_dsxx::SP to an interface type that contains the
 * interface specification defined by rtt_imc::Interface and
 * rtt_imc::Flat_Data_Interface. 
 */
template<class MT>
template<class IT>
Flat_Mat_State_Builder<MT>::Flat_Mat_State_Builder(rtt_dsxx::SP<IT> interface)
    : Mat_State_Builder<MT>()
{
    Require (interface);

    // assign data members from the interface parser
    density            = interface->get_density();
    absorption_opacity = interface->get_absorption_opacity();
    scattering_opacity = interface->get_scattering_opacity();
    temperature        = interface->get_temperature();
    specific_heat      = interface->get_specific_heat();
    implicitness       = interface->get_implicitness_factor();
    delta_t            = interface->get_delta_t();

    Ensure (delta_t > 0.0);
    Ensure (implicitness >= 0.0 && implicitness <= 1.0);
    Ensure (density.size() == absorption_opacity.size());
    Ensure (density.size() == scattering_opacity.size());
    Ensure (density.size() == temperature.size());
    Ensure (density.size() == specific_heat.size());
}

} // end namespace rtt_imc

#endif                          // __imc_Flat_Mat_State_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Flat_Mat_State_Builder.hh
//---------------------------------------------------------------------------//
