//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Comm_Patterns.hh
 * \author Thomas M. Evans
 * \date   Mon May  1 19:50:44 2000
 * \brief  Comm_Patterns header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Comm_Patterns_hh__
#define __mc_Comm_Patterns_hh__

#include "Topology.hh"
#include "ds++/SP.hh"
#include <vector>
#include <map>
#include <utility>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class Comm_Patterns

 * \brief Determine boundary cell communication patterns.

 * The Comm_Patterns class determines the boundary cell relationships between
 * domain discretized, spatial topologies.  In essence, the Comm_Patterns
 * give the following information:
 
 * \arg what processors this processor shares boundary cells with;
 
 * \arg what cells this processor contains that another processor needs
 * information about

 * \arg what cells another processor contains that this processor needs
 * information about

 * Essentially, the Comm_Patterns defines information that the processor
 * needs to do \i gather/scatter type operations.  This type of information
 * is only required in spatially domain decomposed topologies.  The
 * information contained in this class is closely related to the information
 * contained in the Topology classes.  However, the information contained
 * here is related to very specific topology-like information. Also, parallel
 * communication is required to generate the information; thus it belongs in
 * its own class.

 * Accessor functions to Comm_Patterns are in the form of constant iterators
 * to std::map data.  As such, the user can access the data by iterating
 * through it.  The std::map::iterator is a std::pair so the processor is
 * accessed by iterator->first and the cell list by iterator->second[i].
 * Because all accessor functions to this data in Comm_Patterns are const,
 * only map::const_iterator types are returned.  This ensures that clients
 * cannot modify the Comm_Patterns data through the iterators.  However, the
 * iterators themselves can be modified (ie Comm_Patterns::const_iterator i;
 * i++; is valid).  In summary, all accessing of the Comm_Patterns must be
 * done through const_iterator; expressions of the type
 * Comm_Patterns::iterator i = object.get_recv_begin() will be flagged as an
 * error at compile time.  See the examples for usage.

 */
/*!
 * \example mc/test/tstComm_Patterns.cc
 
 * Example usage of Comm_Patterns class.  The full DD mesh is a nine cell
 * mesh (3 x 3 with structured ordering) with cells 1 and 2 on processor 0;
 * cells 3 and 4 on processor 1; cells 5 and 6 on processor 2; and cells 7,
 * 8, and 9 on processor 4.

 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Comm_Patterns 
{
  public:
    // Typedefs.
    typedef std::vector<int>         sf_int;
    typedef std::map<int, sf_int>    proc_map;
    typedef std::pair<int, sf_int>   proc_pair;
    typedef proc_map::iterator       iterator;
    typedef proc_map::const_iterator const_iterator;
    typedef rtt_dsxx::SP<Topology>   SP_Topology;

  private:
    // Map of processor/boundary cell data that the processor sends out to
    // other processors.
    proc_map send_query_map;

    // Map of processor/boundary cell data that the processor receives from
    // other processors.
    proc_map recv_query_map;

  private:
    // Implementation
    void calc_recv_map_DD(SP_Topology);
    void calc_send_map_DD(SP_Topology);

    // Empty the query maps.
    void empty_query_maps();

  public:
    //! Constructor.
    Comm_Patterns() {};
    
    // Calculate communication patterns.
    void calc_patterns(SP_Topology);

    // Query to see of the patterns are set.
    operator bool() const;

    // ACCESSORS
    
    //! Get number of processors from which data is received.
    int get_num_recv_procs() const { return recv_query_map.size(); }

    //! Get number of processors to which data is sent.
    int get_num_send_procs() const { return send_query_map.size(); }

    //! Get iterator to beginning of receive data map.
    const_iterator get_recv_begin() const { return recv_query_map.begin(); } 

    //! Get iterator to end of receive data map.
    const_iterator get_recv_end() const { return recv_query_map.end(); }
};

} // end namespace rtt_mc

#endif                          // __mc_Comm_Patterns_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Comm_Patterns.hh
//---------------------------------------------------------------------------//
