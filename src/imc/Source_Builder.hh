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
#include "mc/Parallel_Data_Operator.hh"
#include "rng/Random.hh"
#include "c4/global.hh"
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
 * processor.  The derived classes of Source_Builder are: 
 *
 * \arg rtt_imc::Rep_Source_Builder for full replication sources; \arg
 * rtt_imc::DD_Source_Builder for full domain decomposition sources; \arg
 * rtt_imc::Serial_Source_Builder for serial sources.
 *
 * Rep_Source_Builder defaults to serial when run on one processor.
 * Serial_Source_Builder is simply an optimized version of Rep_Source_Builder
 * for 1 processor (non-communication) compilations.  Thus, we have full
 * replication, full DD, DD/replication, and serial derived classes.
 * Functionality common to all types of source builders exists in the base
 * class as non-virtual functions.
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
 * edits.  Access functions are provided for this data as regular and virtual
 * functions that are defined in the derived classes 
 *
 * All protected data in the Source_Builder base class is local to a
 * processor.  Global quantities are calculated by using the
 * Parallel_Data_Operations::local_to_global mapper and specifying whether
 * the data to be mapped is Data_Replicated or Data_Decomposed.
 * Specifically, when every replicated cell has exactly the same value, the
 * data is replicated.  When replicated cells have a portion of the global
 * value, the data is said to be decomposed.  If no cells are replicated,
 * i.e., the topology parallel scheme is full domain decomposition, all data
 * is distributed (Data_Distibuted).  In the case of full replication,
 * deterministic quantities such as volume emission energy and surface source
 * energy are data-replicated, and stochastic quantities such as the
 * surviving census particle energies are data-decomposed.
 *
 * The following relationships apply to topology and data decomposition for
 * locally calculated data in "replication" topologies: \arg local data +
 * Data_Replication = global data \arg local data + Data_Decomposed = local
 * data
 *
 * Data_Distributed does not apply to "replication".  In full "DD" all data
 * is Data_Distributed, thus the local data on each processor is inherently
 * global for that cell.  A caveat is that global data fields in "DD"
 * topologies may not contain all the global data elements on a given
 * processor.  It may only contain the data that lives on that processor,
 * ie. the data that is mapped to local cells.
 *
 * Random number stream ID fields are set using the rtt_rng::rn_stream
 * variable.  This variable is \b defined to be global in \b all topologies.
 * For some topologies ("DD") this requires communication after random number
 * IDs are assigned to set and to insure that the rtt_rng::rn_stream value is
 * the same on all processors. 
 *
 * Accessors to data that live in the Source_Builder base class returns the
 * values that have been calculated on processor.  It is up to the client to
 * do any topology-dependent mechanics required to get global (problem-wide)
 * data.  The data could be global or local depending upon the topology.
 * Virtual functions are used to access source data that live in derived
 * classes.  Accordingly, this data is also global or local depending upon
 * the topology.  
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

    //! Elapsed time that problem has been run.
    double elapsed_t;

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

    //! Descriptor for each surface source.
    sf_string ss_desc;

    // BASE CLASS IMPLEMENTATION MEMBER FUNCTIONS

    // Calculate volume emission source energies.
    void calc_evol(const Mat_State<MT> &, const Opacity<MT> &);

    // Calculate surface source energies.
    void calc_ess();

  protected:
    // DATA USED BY ALL SOURCE BUILDERS

    // Interface data.

    //! Requested number of source particles.
    int npnom;

    //! Maximum number of source particles.
    int npmax;

    //! Rate of change of number of source particles per shake.
    double dnpdt;

    //! Problem cycle.
    int cycle;

    //! Timestep.
    double delta_t;

    //! Surface source angular distribution.
    std_string ss_dist;
    
    //! Local (on-processor) census.
    SP_Census census;

    // Data through constructor arguments.

    //! Problem topology.
    SP_Topology topology;

    //! Parallel_Data_Operator for parallel data operations.
    rtt_mc::Parallel_Data_Operator parallel_data_op;

    // Source class data.

    //! Number of particles that we want for this cycle.
    int npwant;

    //! Local field of census energies.
    ccsf_double ecen;

    //! Local (initial) census total energy.
    double ecentot;

    //! Local field of census energy weights.
    ccsf_double ew_cen;

    //! Local field of volume emission energies.
    ccsf_double evol;

    //! Local field of net volume emission energies.
    ccsf_double evol_net;
    
    //! Local field of external material volume source energies.
    ccsf_double mat_vol_src;

    //! Local volume emission total energy.
    double evoltot;

    //! Local field of volume emission energy weights.
    ccsf_double ew_vol;

    //! Local external material volume source total energy.
    double mat_vol_srctot;

    //! Local field of surface source energies.
    ccsf_double ess;

    //! Local field of surface source faces.
    ccsf_int ss_face_in_cell;

    //! Local surface source total energy.
    double esstot;

    //! Local field of surface source energy weights.
    ccsf_double ew_ss;

    //! Local starting random number stream ID field for volume emission.
    ccsf_int volrn;

    //! Local starting random number stream ID field for surface source.
    ccsf_int ssrn;

    // IMPLEMENTATION INHERITANCE

    // Calculate source energies for volume emission and surface source.
    void calc_source_energies(const Mat_State<MT> &, const Opacity<MT> &);

    // Calculate the initial census energy.
    void calc_initial_ecen();

    // Calculate the number of source particles.
    void calc_num_src_particles(const double, const ccsf_double &, 
			       ccsf_int &, int &); 

    // Write the local census on each processor.
    void write_initial_census(SP_Mesh, SP_Rnd_Control, const ccsf_int &,
			      const int &, const ccsf_int &);

    // Comb local census.
    void comb_census(SP_Rnd_Control, int &, double &);

  public:
    // Constructor.
    template<class IT>
    Source_Builder(dsxx::SP<IT>, SP_Mesh, SP_Topology);

    // virtual destructor for correct behavior in inheritance chain
    virtual ~Source_Builder() {/*...*/}

    //! Build dsxx::SP to source.
    virtual SP_Source build_Source(SP_Mesh, SP_Mat_State, SP_Opacity,
				   SP_Rnd_Control) = 0;

    //! Build an initial census from initial radiation temperature.
    virtual void calc_initial_census(SP_Mesh, SP_Mat_State, SP_Opacity,
				     SP_Rnd_Control) = 0;

    //! Get the number of cells accessed by the Source_Builder on processor.
    int num_cells() const { return topology->num_cells(C4::node()); }

    // ACCESSOR FOR CENSUS

    //! Get a dsxx::SP to the census.
    SP_Census get_census() const { return census; }

    // ACCESSOR SERVICES TO GET SOURCE INITIALIZATION DATA

    // These functions are primarily for edits and the problem variables on
    // processor.  The client will be required to use the data in an
    // appropriate fashion based on the problem topology, ie. data could be
    // either global or local depending upon the topology.

    //! Get total volume emission energy on processor.
    double get_evoltot() const { return evoltot; } 

    //! Get global material volume emission source energy/cell on processor.
    double get_mat_vol_src(int cell) const { return mat_vol_src(cell); } 

    //! Get total material volume emission source energy on processor.
    double get_mat_vol_srctot() const { return mat_vol_srctot; }

    //! Get net volume emission energy per cell on processor.
    double get_evol_net(int cell) const { return evol_net(cell); }

    //! Get total surface source energy on processor.
    double get_esstot() const { return esstot; }

    //! Get total initial census energy on processor.
    double get_initial_census_energy() const { return ecentot; }

    // Topology dependent data.  This data could be either global or local
    // depending upon the topology.  The data that is returned lives in the
    // derived class.

    //! Get energy loss in volume emission - topology dependent.
    virtual double get_eloss_vol() const = 0;

    //! Get total number of volume emission particles - toplogy dependent.
    virtual int get_nvoltot() const = 0;

    //! Get number of volume emission particles / cell - topology dependent.
    virtual int get_nvol(int) const = 0;
    
    //! Get energy loss in surface source - topology dependent.
    virtual double get_eloss_ss() const = 0;

    //! Get total number of surface source particles - topology dependent.
    virtual int get_nsstot() const = 0;

    //! Get energy loss in census - topology dependent.
    virtual double get_eloss_cen() const = 0;
    
    //! Get total number of post-comb census particles - topology dependent.
    virtual int get_ncentot() const = 0;
};

} // end namespace rtt_imc

#endif                          // __imc_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Builder.hh
//---------------------------------------------------------------------------//
