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
 * Source_Builder use the rtt_mc::Topology classes to determine where things
 * are and where they need to go.
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
    typedef std::vector<int>                               sf_int;
    typedef std::vector<double>                            sf_double;
    typedef std::vector<std::string>                       sf_string;
    typedef std::string;                                   std_string;
    typedef typename MT::CCSF_double                       ccsf_double;
    typedef typename MT::CCSF_int                          ccsf_int;
    typedef typename MT::CCVF_double                       ccvf_double;

  protected:
    // DATA USED BY ALL SOURCE BUILDERS

    //! Local (on-processor) census.
    SP_Census census;

    //! Local field of census energies.
    ccsf_double ecen;

    //! Local field of census numbers.
    ccsf_int ncen;

  public:
    // Constructor.
    Source_Builder() {}

    // virtual destructor for correct behavior in inheritance chain
    virtual ~Source_Builder() {/*...*/}

    //! Build dsxx::SP to source.
    virtual SP_Source build_Source() = 0;

    // Communication function for DBC checks.
    bool check_global_equiv(int) const;
    bool check_global_equiv(double, double = 1.0e-8) const;

    // ACCESSOR SERVICES TO GET SOURCE INITIALIZATION DATA

    //! Get a dsxx::SP to the census.
    SP_Census get_census() const { return census; }
};

} // end namespace rtt_imc

#endif                          // __imc_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Builder.hh
//---------------------------------------------------------------------------//
