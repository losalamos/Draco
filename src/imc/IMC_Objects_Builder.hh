//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/IMC_Objects_Builder.hh
 * \author Thomas M. Evans
 * \date   Fri Aug 22 16:44:10 2003
 * \brief  IMC_Objects_Builder class definition.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_imc_IMC_Objects_Builder_hh
#define rtt_imc_IMC_Objects_Builder_hh

#include "ds++/SP.hh"
#include "Mat_State.hh"
#include "Opacity.hh"
#include "Diffusion_Opacity.hh"
#include "Source.hh"
#include "Source_Builder.hh"
#include "Tally.hh"
#include "Random_Walk.hh"
#include "Extrinsic_Surface_Tracker.hh"

// Forward declarations.

namespace rtt_rng
{

class Rnd_Control;

}

namespace rtt_mc
{

class Topology;
class Comm_Patterns;

}

namespace rtt_imc
{

class Flat_Data_Interface;
class CDI_Data_Interface;
template<class MT, class FT> class Mat_State_Builder;

//===========================================================================//
/*!
 * \class IMC_Objects_Builder
 * \brief Provides protected functions that build IMC objects during in-cycle
 * initialization.
 *
 * This class provides functionality to build the following, on-processor,
 * IMC objects that are needed to do a cycle of IMC transport:
 * - rtt_imc::Gray_Frequency or rtt_imc::Multigroup_Frequency (FT)
 * - rtt_imc::Mat_State
 * - rtt_imc::Opacity
 * - rtt_imc::Diffusion_Opacity (if needed)
 * - rtt_imc::Source
 * - rtt_imc::Source_Builder (for initialization and edit steps)
 * - rtt_imc::Tally (with appropriate sub tallies)
 * - rtt_imc::Extrinsic_Surface_Tracker (if needed)
 * - rtt_imc::Random_Walk (if needed)
 * .
 * These objects constitute the set of physics objects needed by
 * rtt_imc::Transporter to perform IMC transport for a timestep.  All of
 * these objects are built on-processor; in other words, the interface
 * provides the appropriate data to each on-processor object.
 *
 * The template parameters are Interface Type (IT), Mesh Type (MT), Frequency
 * Type (FT), and Particle Type (PT).  The Interface Type must derive from
 * the following virtual interfaces:
 * - rtt_imc::Interface
 * - rtt_imc::Surface_Tracking_Interface
 * - rtt_imc::CDI_Data_Interface or rtt_imc::Flat_Data_Interface
 * . 
 * In the last case, the interface must only derived from one of the data
 * interfaces.  It parses on the appropriate data builder based on the
 * interface.
 *
 * The services in this class are provided as protected members that can be
 * used by a client to build the appropriate series of IMC objects in a
 * consistent manner.
 *
 * \sa IMC_Objects_Builder.t.hh for detailed descriptions.
 *
 * Code Sample:
 * \code
 *     IMC_Objects_Builder<MT,FT,PT> builder(interface);
 *     builder.build_IMC_objects(mesh, topology, comm_patterns, rnd_control);
 *     frequency = builder.get_Frequency();
 *     mat_state = builder.get_Mat_State();
 *     // ....
 * \endcode
 */
/*! 
 * \example imc/test/tstIMC_Objects_Builder.cc
 * 
 * Test of IMC_Objects_Builder.
 */
// revision history:
// -----------------
// 0) (Fri Aug 22 16:44:10 2003) Thomas M. Evans: original
// 
//===========================================================================//

template<class IT, class MT, class FT, class PT>
class IMC_Objects_Builder 
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<IT>                        SP_Interface;
    typedef rtt_dsxx::SP<MT>                        SP_Mesh;
    typedef rtt_dsxx::SP<FT>                        SP_Frequency;
    typedef rtt_dsxx::SP<Mat_State<MT> >            SP_Mat_State;
    typedef rtt_dsxx::SP<Opacity<MT,FT> >           SP_Opacity;
    typedef rtt_dsxx::SP<Diffusion_Opacity<MT> >    SP_Diff_Opacity;
    typedef rtt_dsxx::SP<Source<MT,FT,PT> >         SP_Source;
    typedef rtt_dsxx::SP<Tally<MT> >                SP_Tally;
    typedef rtt_dsxx::SP<Random_Walk<MT> >          SP_Random_Walk;
    typedef rtt_dsxx::SP<Extrinsic_Surface_Tracker> SP_Tracker;
    typedef rtt_dsxx::SP<rtt_mc::Topology>          SP_Topology;
    typedef rtt_dsxx::SP<rtt_mc::Comm_Patterns>     SP_Comm_Patterns;
    typedef rtt_dsxx::SP<rtt_rng::Rnd_Control>      SP_Rnd_Control;
    typedef rtt_dsxx::SP<Source_Builder<MT,FT,PT> > SP_Source_Builder; 
    typedef rtt_dsxx::SP<Mat_State_Builder<MT,FT> > SP_Mat_State_Builder;

  private: 
    // >>> IMPLEMENTATION

    // Build a CDI_Mat_State_Builder from a CDI_Data_Interface.
    template<class Stop_Explicit_Instantiation>
    SP_Mat_State_Builder get_mat_builder(const CDI_Data_Interface &);

    // Build a Flat_Mat_State_Builder from a Flat_Data_Interface.
    template<class Stop_Explicit_Instantiation>
    SP_Mat_State_Builder get_mat_builder(const Flat_Data_Interface &);

  protected:
    // >>> DATA

    // Interface.
    SP_Interface interface;

    // Built objects.
    
    // Frequency.
    SP_Frequency frequency;
   
    // Material state.
    SP_Mat_State mat_state;

    // Opacity.
    SP_Opacity opacity;
    
    // Diffusion opacity.
    SP_Diff_Opacity diff_opacity;
    
    // Source.
    SP_Source source;

    // Tally.
    SP_Tally tally;

    // Random walk.
    SP_Random_Walk random_walk;

    // Surface tracker.
    SP_Tracker tracker;

    // Builders needed for updates and edits.

    // Source builder (needed to update census energy fields between
    // timesteps). 
    SP_Source_Builder source_builder;

  protected:
    // >>> DERIVED IMPLEMENTATION

    // Build the material state objects.
    void build_mat_state_objects(SP_Mesh);

    // Build the source objects.
    void build_source_builder_object(SP_Mesh, SP_Topology, SP_Rnd_Control);
    void build_source_object(SP_Mesh, SP_Rnd_Control, SP_Comm_Patterns);

    // Build the tally object.
    void build_tally_object(SP_Mesh);

    // Build particle transport objects.
    void build_time_dependent_particle_objects(SP_Mesh);
    void build_time_independent_particle_objects(SP_Mesh);

    // Reset a member object.
    template<class T>
    void reset(rtt_dsxx::SP<T> &object) { object = rtt_dsxx::SP<T>(); }

  public:
    // Constructor.
    explicit IMC_Objects_Builder(SP_Interface);

    // Destructor.
    virtual ~IMC_Objects_Builder() {/*...*/}

    // >>> PUBLIC INTERFACE

    // >>> ACCESSORS

    //! Get the FT (rtt_imc::Gray_Frequency or rtt_imc::Multigroup_Frequency).
    SP_Frequency get_Frequency() const { return frequency; }

    //! Get the material state (rtt_imc::Mat_State).
    SP_Mat_State get_Mat_State() const { return mat_state; }

    //! Get the opacity (rtt_imc::Opacity).
    SP_Opacity get_Opacity() const { return opacity; }
 
    //! Get the diffusion opacity (rtt_imc::Diffusion_Opacity). 
    SP_Diff_Opacity get_Diffusion_Opacity() const { return diff_opacity; }

    //! Get the source (rtt_imc::Source).
    SP_Source get_Source() const { return source; }

    //! Get the tally (rtt_imc::Tally).
    SP_Tally get_Tally() const { return tally; }

    //! Get the random walk object (rtt_imc::Random_Walk).
    SP_Random_Walk get_Random_Walk() const { return random_walk; }

    //! Get the surface tracker (rtt_imc::Extrinsic_Surface_Tracker).
    SP_Tracker get_Surface_Tracker() const { return tracker; }

    //! Get the source builder (rtt_imc::Source_Builder).
    SP_Source_Builder get_Source_Builder() const { return source_builder; }
};

} // end namespace rtt_imc

#endif // rtt_imc_IMC_Objects_Builder_hh

//---------------------------------------------------------------------------//
//              end of imc/IMC_Objects_Builder.hh
//---------------------------------------------------------------------------//
