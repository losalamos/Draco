//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Rep_Source_Builder.hh
 * \author Thomas M. Evans
 * \date   Thu Dec  9 10:31:16 1999
 * \brief  Header file for Rep_Source_Builder.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Rep_Source_Builder_hh__
#define __imc_Rep_Source_Builder_hh__

#include "Source_Builder.hh"
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
 * \class Rep_Source_Builder
 *
 * \brief Rep_Source_Builder is a derived rtt_imc::Source_Builder class for
 * full replication problems.
 *
 * Rep_Source_Builder is used to build rtt_imc::Source objects in full
 * replication topologies.
 *
 * As with the Source_Builder base class, Rep_Source_Builder is templated on
 * Mesh Type (MT) and Particle Type (PT).  The PT defaults to
 * rtt_imc::Particle.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT, class PT = Particle<MT> >
class Rep_Source_Builder : public Source_Builder<MT,PT>
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

  private:
    //! Build an initial census.
    void calc_initial_census();

  public:
    // Constructor.
    Rep_Source_Builder();

    // Build source.
    SP_Source build_Source();
};

} // end namespace rtt_imc

#endif                          // __imc_Rep_Source_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Rep_Source_Builder.hh
//---------------------------------------------------------------------------//
