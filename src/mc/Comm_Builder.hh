//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Comm_Builder.hh
 * \author Thomas M. Evans
 * \date   Tue Jun  1 17:01:44 1999
 * \brief  Comm_Builder header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Comm_Builder_hh__
#define __mc_Comm_Builder_hh__

#include "Communicator.hh"
#include "mc/Topology.hh"
#include "ds++/SP.hh"
#include <vector>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class Comm_Builder

 * \brief Build a rtt_mc::Communicator using the rtt_mc::Topology class.

 * The Comm_Builder builds rtt_mc::Communicator objects using the
 * rtt_mc::Topology class.  Because the Topology class contains all the
 * information needed by the Communicator to determine where particles go
 * when they cross processor boundaries, the Comm_Builder does not have to
 * perform any communication.  Thus, the Communicator is built on each
 * processor in parallel.
 
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class PT>
class Comm_Builder 
{
  public:
    // Usefull typedefs in Comm_Builder.
    typedef rtt_dsxx::SP<Communicator<PT> > SP_Communicator;
    typedef rtt_dsxx::SP<rtt_mc::Topology>  SP_Topology;
    typedef std::vector<int>                sf_int;
    typedef std::vector<std::vector<int> >  vf_int;

  private:
    // Build a communicator in full DD topologies.
    SP_Communicator build_DD_Comm(SP_Topology);

  public:
    // default constructor
    Comm_Builder() {}

    // Deprecated build service for use with certain host codes. 
    SP_Communicator build_Communicator(const sf_int &, const sf_int &);

    // Standard build service.
    SP_Communicator build_Communicator(SP_Topology);
};

} // end namespace rtt_mc

#endif                          // __mc_Comm_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Comm_Builder.hh
//---------------------------------------------------------------------------//
