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
 * The Mat_State_Builder provides a simple interface to build two objects:
 *
 * \arg \b build_Mat_State() builds a Mat_State object
 * \arg \b build_Frequency() builds a Frequency object
 * \arg \b build_Opacity() builds an Opacity object
 * \arg \b build_Frequency_Operations() builds a Frequency_Operations object
 *
 * The Mat_State must be built first because the Mat_State object is an
 * argument to build_Opacity().  The argument is required for efficiency's
 * sake; we do not want to reread specific heats, etc. from an input file to
 * build the opacities after we have already gotten them to build the
 * Mat_State.
 *
 * See the tstMat_State_Builder.cc for usage examples.
 */
/*!
 * \example imc/test/tstMat_State_Builder.cc
 *
 * Test of Mat_State_Builder types and rtt_imc::Opacity,
 * rtt_imc::Diffusion_Opacity, and rtt_imc::Mat_State.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class FT>
class Mat_State_Builder 
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<MT>               SP_Mesh;
    typedef rtt_dsxx::SP<Mat_State<MT> >   SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >  SP_Opacity;
    typedef rtt_dsxx::SP<FT>               SP_Frequency;
    
  public:
    //! Constructor.
    Mat_State_Builder() {/*...*/}

    //! Virtual destructor.
    virtual ~Mat_State_Builder() {/* need a destructor for inheritance */}

    // >>> PUBLIC INTERFACE

    //! Build a frequency definition.
    virtual SP_Frequency build_Frequency() = 0;

    //! Build a material state.
    virtual SP_Mat_State build_Mat_State(SP_Mesh) = 0;

    //! Build an opacity.
    virtual SP_Opacity build_Opacity(SP_Mesh, SP_Frequency, SP_Mat_State) = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Mat_State_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Mat_State_Builder.hh
//---------------------------------------------------------------------------//
