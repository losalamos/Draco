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
 * processor.  Thus, we have full replication, full DD, DD/replication, and
 * serial derived classes.  Functionality common to all types of source
 * builders exists in the base class as non-virtual functions.
 *
 * Source_Builder and Source_Builder hierarchies are templated on mesh and
 * particle type.  By default, the particle type is rtt_imc::Particle.
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
    // Constructor.
    Source_Builder() {}

    // virtual destructor for correct behavior in inheritance chain
    // virtual ~Source_Builder() {/*...*/}

    // Communication function for DBC checks.
    bool check_global_equiv(int) const;
    bool check_global_equiv(double, double = 1.0e-8) const;
};

} // end namespace rtt_imc

#endif                          // __imc_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Source_Builder.hh
//---------------------------------------------------------------------------//
