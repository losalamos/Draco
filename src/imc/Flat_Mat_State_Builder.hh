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
#include "Flat_Data_Container.hh"
#include "Mat_State.hh"
#include "Opacity.hh"
#include "Diffusion_Opacity.hh"
#include "Frequency.hh"
#include "Global.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
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

template<class MT, class FT>
class Flat_Mat_State_Builder : public Mat_State_Builder<MT,FT>
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>                                   SP_Mesh;
    typedef rtt_dsxx::SP<Mat_State<MT> >                       SP_Mat_State;
    typedef rtt_dsxx::SP<FT>                                   SP_Frequency;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >                      SP_Opacity;
    typedef rtt_dsxx::SP<Diffusion_Opacity<MT> >               SP_Diff_Opacity;
    typedef std::vector<double>                                sf_double;
    typedef std::vector<sf_double>                             vf_double;
    typedef rtt_imc::global::Type_Switch<Gray_Frequency>       Switch_Gray;
    typedef rtt_imc::global::Type_Switch<Multigroup_Frequency> Switch_MG;
    typedef rtt_dsxx::SP<Flat_Data_Container>                  SP_Flat_Data;
    typedef rtt_dsxx::SP<Opacity<MT,Gray_Frequency> >          SP_Gray_Opacity;
    typedef rtt_dsxx::SP<Opacity<MT,Multigroup_Frequency> >    SP_MG_Opacity;
    typedef rtt_dsxx::SP<Gray_Frequency>                       SP_Gray;
    typedef rtt_dsxx::SP<Multigroup_Frequency>                 SP_MG;
    typedef typename rtt_imc::global::Type_Switch<FT>::Type    Dummy_Type;

  private:
    // Flat, cell-centered data fields received from the interface.
    SP_Flat_Data flat_data;

    // Densities in g/cc.
    sf_double    density;

    // Material temperatures in keV.
    sf_double    temperature;

    // Fleck and Cummings implicitness factor.
    double       implicitness;

    // Timestep in shakes.
    double       delta_t;

  private:
    // >>> BUILT OBJECTS
    
    // Built frequency.
    SP_Frequency frequency;

    // Built Mat_State.
    SP_Mat_State mat_state;

    // Built Opacity.
    SP_Opacity opacity;
    
    // Built Diffusion_Opacity.
    SP_Diff_Opacity diff_opacity;

  private:
    // >>> IMPLEMENTATION

    // Build the Mat_State.
    void build_Mat_State(SP_Mesh);

  private:
    // >>> PARTIAL SPECIALIZATIONS ON FREQUENCY TYPE

    // Build a Gray_Frequency.
    template<class Stop_Explicit_Instantiation>
    rtt_dsxx::SP<Gray_Frequency> build_frequency(Switch_Gray);

    // Build a Multigroup_Frequency
    template<class Stop_Explicit_Instantiation>
    rtt_dsxx::SP<Multigroup_Frequency> build_frequency(Switch_MG);

    // Build an Opacity<MT,Gray_Frequency>
    template<class Stop_Explicit_Instantiation>
    SP_Gray_Opacity build_opacity(Switch_Gray, SP_Mesh);

    // Build an Opacity<MT,Multigroup_Frequency>
    template<class Stop_Explicit_Instantiation>
    SP_MG_Opacity build_opacity(Switch_MG, SP_Mesh);

  public:
    // Constructor.
    template<class IT>
    explicit Flat_Mat_State_Builder(rtt_dsxx::SP<IT>);

    // >>> PUBLIC INTERFACE

    // Build frequency, material state, and opacity.
    void build_mat_classes(SP_Mesh);

    //! Get the frequency.
    SP_Frequency get_Frequency() const { return frequency; }
    
    //! Get the mat state.
    SP_Mat_State get_Mat_State() const { return mat_state; }

    //! Get the opacity
    SP_Opacity get_Opacity() const { return opacity; }

    //! Get the diffusion opacity.
    SP_Diff_Opacity get_Diffusion_Opacity() const { return diff_opacity; }
};

//---------------------------------------------------------------------------//
// MEMBER TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 *
 * The constructor gets data from the interface type.  The interface type
 * must provide member functions defined by rtt_imc::Interface and
 * rtt_imc::Flat_Data_Interface.  However, because this is used as a template
 * argument, it is not required that the interface type inherit from these
 * classes.  It is only required that the interface type have these functions
 * defined.
 *
 * \param interface rtt_dsxx::SP to an interface type that contains the
 * interface specification defined by rtt_imc::Interface and
 * rtt_imc::Flat_Data_Interface. 
 */
template<class MT, class FT>
template<class IT>
Flat_Mat_State_Builder<MT,FT>::Flat_Mat_State_Builder(
    rtt_dsxx::SP<IT> interface)
    : Mat_State_Builder<MT,FT>()
{
    Require (interface);

    // assign data members from the interface parser
    density       = interface->get_density();
    temperature   = interface->get_temperature();
    implicitness  = interface->get_implicitness_factor();
    delta_t       = interface->get_delta_t();
    flat_data     = interface->get_flat_data_container();

    Ensure (delta_t > 0.0);
    Ensure (implicitness >= 0.0 && implicitness <= 1.0);
    Ensure (density.size() == temperature.size());
    Ensure (flat_data);
    Ensure (flat_data->specific_heat.size() == density.size());
}

} // end namespace rtt_imc

#endif                          // __imc_Flat_Mat_State_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Flat_Mat_State_Builder.hh
//---------------------------------------------------------------------------//
