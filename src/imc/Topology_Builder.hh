//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Topology_Builder.hh
 * \author Thomas M. Evans
 * \date   Mon Nov 29 14:33:42 1999
 * \brief  Topology_Builder header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_Topology_Builder_hh__
#define __imc_Topology_Builder_hh__

#include "mc/Topology.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_imc
{
 
//===========================================================================//
/*!
 * \class Topology_Builder
 *
 * \brief Topology_Builder builds a rtt_mc::Topology object.
 *
 * The Topology_Builder class is used to build rtt_mc::Topology objects based
 * on various spatial mesh decompositions.  As such, the Topology_Builder is
 * templated on the mesh type.  The Topology_Builder constructor requires an
 * instantiation of the global mesh object to determine the complete parallel
 * topology for the problem.  This mesh object may be destroyed after usage.
 *
 * At present, the Topology_Builder can do two types of decomposition,
 * replication and DD.  Replication is where the mesh is fully replicated on
 * each processor.  DD is where each mesh cell is uniquely placed on a
 * certain processor.  To perform this topology construction, the
 * Topology_Builder determines the processor capacity from an appropriate
 * interface.  If the capacity is greater than the number of cells then full
 * replication is performed.  Otherwise, full DD is performed.
 *
 * Full DD topologies are instantiated as rtt_mc::Topology smart pointers to
 * a rtt_mc::General_Topology instance.  Replication topologies are
 * instantiated as rtt_mc::Topology smart pointers to a rtt_mc::Rep_Topology
 * instance.
 *
 * In the future, the Topology_Builder will support more advanced
 * (complicated) parallel decompositions.  In these cases, the
 * Topology_Builder could take source initialization information to partition
 * the mesh.
 *
 * The Topology_Builder is designed to work on the host processor during
 * problem setup.  However, in the future, one could easily add alternate
 * constructors and build functions to work in a multi-processor setup mode.
 */
/*!
 * \example imc/test/tstTopology_Builder.cc
 *
 * Example usage of the Topology_Builder class.  In this example Topologies
 * are built for both full replication and full DD based on simple interface
 * classes.  The mesh is a six cell 2D mesh.
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

template<class MT>
class Topology_Builder 
{
  public:
    // Usefull typedefs in this class.
    typedef rtt_dsxx::SP<rtt_mc::Topology> SP_Topology;
    typedef rtt_dsxx::SP<MT>               SP_MT;
    typedef std::vector<std::vector<int> > vf_int;
    typedef std::string                    std_string;

  private:
    // The cells/processor capacity.
    int capacity;

    // build a full DD topology
    void build_DD(vf_int &, vf_int &, vf_int &, std_string &, SP_MT);

  public:
    // Constructor, templated on the Interface Type.
    template<class IT> Topology_Builder(rtt_dsxx::SP<IT>);

    // Builder functions.
    SP_Topology build_Topology(SP_MT);
};

} // end namespace rtt_imc

#endif                          // __imc_Topology_Builder_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Topology_Builder.hh
//---------------------------------------------------------------------------//
