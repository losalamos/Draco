//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Source_Builder.hh
 * \author Thomas M. Evans
 * \date   Wed Dec  8 14:35:42 1999
 * \brief  Source_Builder header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Source_Builder_hh__
#define __imc_Source_Builder_hh__

#include "Particle.hh"
#include "Source.hh"
#include "Particle_Buffer.hh"
#include "Opacity.hh"
#include "Mat_State.hh"
#include "mc/Topology.hh"
#include "rng/Random.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>
#include <iostream>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Source_Builder
 *
 * \brief Source_Builder abstract base class for building IMC sources.
 *
 * Source_Builder is an abstract base class used to create source builder
 * inheritance hierarchies.  In short, Source_Builder builds an IMC source on
 * processor.  The derived classes of Source_Builder are: \arg
 * rtt_imc::Rep_Source_Builder for full replication sources; \arg
 * rtt_imc::DD_Source_Builder for full domain decomposition sources; \arg
 * rtt_imc::Serial_Source_Builder for serial sources.  Rep_Source_Builder
 * defaults to serial when run on one processor.  Serial_Source_Builder is
 * simply an optimized version of Rep_Source_Builder for 1 processor
 * (non-communication) compilations.  Thus, we have full replication, full
 * DD, DD/replication, and serial derived classes.  Functionality common to
 * all types of source builders exists in the base class as non-virtual
 * functions.
 *
 * Source_Builder and Source_Builder hierarchies are templated on mesh and
 * particle type.  By default, the particle type is rtt_imc::Particle.
 *
 * Source_Builder classes use the rtt_mc::Topology classes to determine where
 * things are and where they need to go.
 *
 * Source_Builder is used to build smart pointers to rtt_imc::Source on each
 * processor.  Additionally, the Source_Builder class, in the process of
 * making sources, calculates global source information that can be used for
 * edits.  Access functions are provided for this data.  
 */
/*!
 * \example imc/test/tstSource_Builder.cc
 *
 * Example usage of the Source_Builder class.  The derived classes
 * Rep_Source_Builder and DD_Source_Builder are tested.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class PT = Particle<MT> >
class Source_Builder 
{
  public:
    // typedefs used in the inheritance chain
    typedef dsxx::SP<Source<MT,PT> >                       SP_Source;
    typedef dsxx::SP<typename Particle_Buffer<PT>::Census> SP_Census;
    typedef dsxx::SP<Opacity<MT> >                         SP_Opacity;
    typedef dsxx::SP<Mat_State<MT> >                       SP_Mat_State;
    typedef dsxx::SP<MT>                                   SP_Mesh;
    typedef dsxx::SP<rtt_rng::Rnd_Control>                 SP_Rnd_Control;
    typedef dsxx::SP<rtt_mc::Topology>                     SP_Topology;
    typedef std::vector<int>                               sf_int;
    typedef std::vector<std::vector<int> >                 vf_int;
    typedef std::vector<double>                            sf_double;
    typedef std::vector<std::string>                       sf_string;
    typedef std::string                                    std_string;
    typedef typename MT::CCSF_double                       ccsf_double;
    typedef typename MT::CCSF_int                          ccsf_int;
    typedef typename MT::CCVF_double                       ccvf_double;

  private:
    // BASE CLASS DATA

    // Interface data.

    //! Vector of local external volume source energies.
    sf_double evol_ext;

    //! External radiation source duration.
    double rad_s_tend;

    //! External radiation source energies.
    sf_double rad_source;

    //! Initial radiation temperature, which determine initial census.
    sf_double rad_temp;

    //! Position (and number) of surface sources.
    sf_string ss_pos;

    //! Temperatures of surface sources.
    sf_double ss_temp;

    //! List of local cells that have surface sources.
    vf_int defined_surcells;

    // BASE CLASS IMPLEMENTATION MEMBER FUNCTIONS

    // Calculate volume emission source energies.
    void calc_evol(const Mat_State<MT> &, const Opacity<MT> &);

    // Calculate surface source energies.
    void calc_ess();

  protected:
    // DATA USED BY ALL SOURCE BUILDERS

    // Interface data.

    //! Problem cycle.
    int cycle;

    //! Timestep.
    double delta_t;

    //! Problem topology.
    SP_Topology topology;
    
    //! Local (on-processor) census.
    SP_Census census;

    // Source class data.

    //! Local field of census energies.
    ccsf_double ecen;

    //! Local (initial) census total energy.
    double ecentot;

    //! Local field of census numbers.
    ccsf_int ncen;

    //! Local field of volume emission energies.
    ccsf_double evol;

    //! Local field of net volume emission energies.
    ccsf_double evol_net;
    
    //! Local field of external material volume source energies.
    ccsf_double mat_vol_src;

    //! Local volume emission total energy.
    double evoltot;

    //! Local external material volume source total energy.
    double mat_vol_srctot;

    //! Local field of surface source energies.
    ccsf_double ess;

    //! Local field of surface source faces.
    ccsf_int ss_face_in_cell;

    //! Local surface source total energy.
    double esstot;

    //! Local starting random number stream ID field for volume emission.
    ccsf_int volrn;

    //! Local starting random number stream ID field for surface source.
    ccsf_int ssrn;

    //! Local starting random number stream ID vector for initial census.
    // We use a vector field here because the cenrn field is only required
    // during census initialization.  Thus, after the census has been created 
    // we save memory be setting this variable to zero.
    sf_int cenrn;

    // IMPLEMENTATION INHERITANCE

    // Calculate source energies for volume emission and surface source.
    void calc_source_energies(const Mat_State<MT> &, const Opacity<MT> &);

    // Calculate the initial census energy.
    void calc_ecen_init();

    // Calculate the number of source particles.
    int calc_num_src_particles(const double, const ccsf_double &, 
			       ccsf_int &); 

  public:
    // Constructor.
    template<class IT>
    Source_Builder(dsxx::SP<IT>, SP_Mesh, SP_Topology);

    // virtual destructor for correct behavior in inheritance chain
    virtual ~Source_Builder() {/*...*/}

    //! Build dsxx::SP to source.
    virtual SP_Source build_Source(SP_Mesh, SP_Mat_State, SP_Opacity,
				   SP_Rnd_Control) = 0; 

    // ACCESSOR SERVICES TO GET SOURCE INITIALIZATION DATA

    //! Get a dsxx::SP to the census.
    SP_Census get_census() const { return census; }
};

} // end namespace rtt_imc

#endif                          // __imc_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Builder.hh
//---------------------------------------------------------------------------//
