//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Flat_Mat_State_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:36:42 2001
 * \brief  Flat_Mat_State_Builder class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_Flat_Mat_State_Builder_hh
#define rtt_imc_Flat_Mat_State_Builder_hh

#include "Mat_State_Builder.hh"
#include "Flat_Data_Container.hh"
#include "Mat_State.hh"
#include "Opacity.hh"
#include "Diffusion_Opacity.hh"
#include "Frequency.hh"
#include "Hybrid_Diffusion.hh"
#include "Global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include <vector>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Flat_Mat_State_Builder
 *
 * This class builds instances of rtt_imc::Opacity,
 * rtt_imc::Diffusion_Opacity, and rtt_imc::Mat_State.  It can be templated
 * on either Gray_Frequency or Multigroup_Frequency through the FT (Frequency
 * Type) parameter. It receives its data from a defined interface in a "flat"
 * form.  In other words, the data fields are already defined.  It receives
 * an Interface Type (IT) in its constructor.  The IT must be a derived class
 * of rtt_imc::Interface and rtt_imc::Flat_Data_Interface.
 *
 * For gray problems the Flat_Data_Container::gray_absorption_opacity and
 * Flat_Data_Container::gray_scattering_opacity fields must be allocated.
 * For multigroup problems the following fields must be allocated:
 * - Flat_Data_Container::mg_absorption_opacity
 * - Flat_Data_Container::mg_scattering_opacity
 * - Flat_Data_Container::group_boundaries
 * . 
 * Both gray and multigroup problems require the
 * Flat_Data_Container::specific_heat field.
 *
 * To build Diffusion_Opacity objects for FT = Gray_Frequency, the
 * rtt_imc::Flat_Data_Container::rosseland_opacity field must be filled. As
 * opposed to the rtt_imc::CDI_Mat_State_Builder, the Flat_Mat_State_Builder
 * does not attempt to assemble Rosseland opacities from the absorption and
 * scattering opacities.  The Rosseland opacities must come through the
 * interface in the Flat_Data_Container.
 *
 * For multigroup, the Rosseland opacities are calculated by integrating over
 * each group as follows:
 * \f[
 * \frac{1}{\sigma_{R}} = \frac{\int\frac{1}{(\sigma_{a}+\sigma_{s})}
 * \frac{\partial B}{\partial T}}{\int\frac{\partial B}{\partial T}}
 * \f]  
 *
 * This class, along with rtt_imc::CDI_Mat_State_Builder, is a derived class
 * of Mat_State_Builder.  It should be used when a client already has
 * cell-centered material defined.
 *
 * \sa rtt_imc::Interface, rtt_imc::Flat_Data_Interface,
 * rtt_imc::Mat_State_Builder.  See tstFlat_Mat_Data_Builder for examples of
 * usage.
 */
/*!
 * \example imc/test/tstFlat_Mat_State_Builder.cc
 *
 * Test of CDI_Mat_State_Builder derived type.
 */
// revision history:
// -----------------
// 0) original
// 1) 20-FEB-2003 : updated to use new Mat_State_Builder interface
// 2) 13-MAR-2003 : updated to build diffusion opacities for gray problems
// 3) 28-JUL-2003 : updated to build diffusion opacities for mg problems
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

    // Switch for building diffusion opacities.
    bool          build_diffusion_opacity;

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
    void build_mat_state(SP_Mesh);

  private:
    // >>> PARTIAL SPECIALIZATIONS ON FREQUENCY TYPE

    // Build a Gray_Frequency.
    template<class Stop_Explicit_Instantiation>
    void build_frequency(Switch_Gray);

    // Build a Multigroup_Frequency
    template<class Stop_Explicit_Instantiation>
    void build_frequency(Switch_MG);

    // Build an Opacity<MT,Gray_Frequency>
    template<class Stop_Explicit_Instantiation>
    void build_opacity(Switch_Gray, SP_Mesh);

    // Build an Opacity<MT,Multigroup_Frequency>
    template<class Stop_Explicit_Instantiation>
    void build_opacity(Switch_MG, SP_Mesh);

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
    : Mat_State_Builder<MT,FT>(),
      density(interface->get_density()),
      temperature(interface->get_temperature()),
      implicitness(interface->get_implicitness_factor()),
      delta_t(interface->get_delta_t()),
      flat_data(interface->get_flat_data_container()),
      build_diffusion_opacity(false)
{
    Require (interface);
    Require (flat_data);
    
    // set switch
    int hybrid = interface->get_hybrid_diffusion_method();
    switch (hybrid)
    {
    case Hybrid_Diffusion::TRANSPORT:
	build_diffusion_opacity = false;
	break;

    case Hybrid_Diffusion::RANDOM_WALK:
    case Hybrid_Diffusion::DDIMC:
	build_diffusion_opacity = true;
	break;

    default:
	throw rtt_dsxx::assertion("Invalid hybrid diffusion scheme.");
	break;
    }

    Ensure (delta_t > 0.0);
    Ensure (implicitness >= 0.0 && implicitness <= 1.0);
    Ensure (density.size() == temperature.size());
    Ensure (flat_data);
    Ensure (flat_data->specific_heat.size() == density.size());
}

} // end namespace rtt_imc

#endif                          // rtt_imc_Flat_Mat_State_Builder_hh

//---------------------------------------------------------------------------//
//                              end of imc/Flat_Mat_State_Builder.hh
//---------------------------------------------------------------------------//
