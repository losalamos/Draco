//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Topology.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 17 14:56:25 1999
 * \brief  Topology class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Topology_hh__
#define __mc_Topology_hh__

#include <string>
#include <vector>
#include <iostream>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class Topology
 *
 * \brief Topology is a parallel-topology abstract base class that maps cells
 * to processors.
 *
 * The Topology keeps a record of processor<-->mesh-cell relationships. Thus,
 * the Topology class can tell a client across-processor information about a
 * global cell.  This information includes: 
 * \arg \b global_cell: gives the global_cell index for a local_cell on a
 * given processor
 * \arg \b local_cell: gives the local_cell index for a global_cell and
 * processor index
 * \arg \b processor: gives the processor for a global_cell index
 * \arg \b cell_list: gives the global cells on a given processor
 * \arg \b processor_list: gives the processors that a global_cell is on
 * \arg \b boundary_cells: returns boundary cell information for a processor.
 *
 * Topology has two derived classes-a rtt_mc::General_Topology class that is
 * used for DD and DD/replication topologies and a rtt_mc::Rep_Topology class
 * that is used for full replication topologies.  The General_Topology class
 * can also be applied to full replication problems; however, because this
 * topology is so simple, this is a waste of memory and effort.  All
 * functionality in the base class is available through public inheritance in
 * the derived classes.
 *
 * The topology family of classes should be build by an appropriate builder
 * class.  The constructors for the derived classes demand different
 * information.  See rtt_mc::General_Topology and rtt_mc::Rep_Topology for
 * the appropriate building requirements. 
 */
/*!
 * \example mc/test/tstTopology.cc
 *
 * Example usage of the Topology class.  In this example the Topology class
 * is built directly.  In general, Topology instances will be built by
 * Topology builders that take data from an interface and perform some
 * particular spatial partitioning.
 */
// revision history:
// -----------------
// 0) original
// 1) 30-NOV-99 : made Topology an inheritance tree --> derived classes are
//                accessed through polymorphism to describe general
//                topologies, replication topologies, and local topologies.
// 
//===========================================================================//

class Topology 
{
  public:
    // typedefs used in this class
    typedef std::vector<int>               sf_int;
    typedef std::vector<std::vector<int> > vf_int;
    typedef std::string                    std_string;

  private:
    // parallel scheme for this topology
    std_string parallel_scheme;

  public:
    // Constructor
    Topology(const std_string &);

    //! Destructor for proper performance in inheritance tree.
    virtual ~Topology() {/*...*/}

    //! Return the parallel scheme.
    const std_string& get_parallel_scheme() const { return parallel_scheme; }
    
    // virtual functions shared by derived classes

    //! Return the number of global cells.
    virtual int num_cells() const = 0;

    //! Return the number of cells on processor.
    virtual int num_cells(int) const = 0;

    //! Return the number of processors storing a global cell.
    virtual int num_procs(int) const = 0;

    //! Get the global cell index on processor.
    virtual int global_cell(int) const = 0;

    //! Get the global cell index on a given processor.
    virtual int global_cell(int, int) const = 0;

    //! Get the local cell index on processor.
    virtual int local_cell(int) const = 0;
    
    //! Get the local cell index on a given processor.
    virtual int local_cell(int, int) const = 0;

    //! Get a list of global cells on processor.
    virtual sf_int get_cells(int) const = 0;

    //! Get a list of processors that a global cell is on.
    virtual sf_int get_procs(int) const = 0;

    //! Diagnostic printing.
    virtual void print(std::ostream &) const = 0;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

//! Overloaded stream output operator.
std::ostream& operator<<(std::ostream &, const Topology &);

} // end namespace rtt_mc

#endif                          // __mc_Topology_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Topology.hh
//---------------------------------------------------------------------------//
