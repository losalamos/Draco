//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Mat_State_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 15:59:23 2001
 * \brief  Mat_State_Builder class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Mat_State_Builder_hh__
#define __imc_Mat_State_Builder_hh__

#include "ds++/SP.hh"

namespace rtt_imc
{

// Forward declarations
template<class MT> class Mat_State;
template<class MT, class FT> class Opacity;
template<class MT> class Diffusion_Opacity;

//===========================================================================//
/*!
 * \class Mat_State_Builder
 *
 * Base class definitions for material data builders.  The material data
 * builders build the Opacity and Mat_State classes each timestep.
 * Currently, there are two Mat_State_Builder derived classes:
 *
 * \arg \b Flat_Mat_State_Builder: accepts data through a flat interface
 * defined by the rtt_imc::Flat_Data_Interface.
 *
 * \arg \b CDI_Mat_State_Builder: uses the rtt_cdi package to read data;
 * necessary input is defined by the rtt_imc::CDI_Data_Interface
 *
 * The rtt_imc::Interface class defines additional data used by the
 * Mat_State_Builder derived classes.
 *
 * The class is templated on mesh type.
 *
 * The Mat_State_Builder provides a simple interface to build the following
 * objects: 
 *
 * \arg \b rtt_imc::Gray_Frequency gray frequency group layout
 * \arg \b rtt_imc::Multigroup_Frequency multigroup frequency group layout
 * \arg \b rtt_imc::Mat_State material state data
 * \arg \b rtt_imc::Opacity opacities used by transport
 * \arg \b rtt_imc::Diffusion_Opacity opacities used by diffusion
 *
 * The objects are built by a call to \p build_mat_classes(). The built
 * objects can be retrieved from the builder by calling
 * - get_Frequency()
 * - get_Mat_State()
 * - get_Opacity()
 * - get_Diffusion_Opacity()
 * .
 *
 * The rtt_imc::Diffusion_Opacity class is built if the
 * rtt_imc::Interface::get_hybrid_diffusion_method function returns an
 * appropriate method requiring diffusion opacities.  These methods are
 * defined in the rtt_imc::Hybrid_Diffusion::Method enumeration.
 *
 * \sa See the tstCDI_Mat_State_Builder.cc and tstFlat_Mat_State_Builder.cc
 * for usage examples.
 */
// revision history:
// -----------------
// 0) original
// 1) 20-FEB-2003 : updated public interface
// 
//===========================================================================//

template<class MT, class FT>
class Mat_State_Builder 
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>                     SP_Mesh;
    typedef rtt_dsxx::SP<Mat_State<MT> >         SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >        SP_Opacity;
    typedef rtt_dsxx::SP<Diffusion_Opacity<MT> > SP_Diffusion_Opacity;
    typedef rtt_dsxx::SP<FT>                     SP_Frequency;
    
  public:
    //! Constructor.
    Mat_State_Builder() {/*...*/}

    //! Virtual destructor.
    virtual ~Mat_State_Builder() {/* need a destructor for inheritance */}

    // >>> BUILD INTERFACE

    //! Build frequency, material state, and opacity.
    virtual void build_mat_classes(SP_Mesh) = 0;

    //! Get the frequency.
    virtual SP_Frequency get_Frequency() const = 0;
    
    //! Get the mat state.
    virtual SP_Mat_State get_Mat_State() const = 0;

    //! Get the opacity
    virtual SP_Opacity get_Opacity() const = 0;

    //! Get the diffusion opacity.
    virtual SP_Diffusion_Opacity get_Diffusion_Opacity() const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Mat_State_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Mat_State_Builder.hh
//---------------------------------------------------------------------------//
