//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/General_Topology.hh
 * \author Thomas M. Evans
 * \date   Tue Nov 30 13:19:05 1999
 * \brief  General_Topology header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_General_Topology_hh__
#define __mc_General_Topology_hh__

#include "Topology.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <algorithm>

namespace rtt_mc
{
 
//===========================================================================//
/*!
 * \class General_Topology
 *
 * \brief General_Topology is a derived rtt_mc::Topology class for general
 * spatial topologies.
 *
 * General_Topology is a fully general Topology class that handles all
 * topology types: "replication", "DD", and "DD/replication".  However, it is
 * bulky for "replication" topologies.  Fully replicated problem domains
 * should use rtt_mc::Rep_Topology.
 *
 * General_Topology needs to be built by a Topology_Builder that may vary
 * from code to code depending upon how the topology is constructed.  The
 * General_Topology constructor takes complete spatial state information as
 * arguments.  This information is constructed based upon the problem setup,
 * source, mesh, etc.  The General_Topology class constructor requires three
 * vector<vector<int>> fields that contain the following information:
 * \arg \b cells_per_proc: global cells per processor.
 * \arg \b procs_per_cell: processors per global cell.
 * \arg \b bound_cells: list of boundary cells per processor.
 *
 * Additionally, a string is provided that states the parallel topology
 * description that is one of the following:
 * \arg \b "replication": full replication of cells across processors
 * \arg \b "DD": full domain decomposition of cells across processors
 * \arg \b "DD/replication": hybrid scheme where some cells are replicated 
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class General_Topology : public Topology
{
  private:
    // local_cells per processor list
    vf_int cells_per_proc;
    
    // processors per global_cell list
    vf_int procs_per_cell;

    // boundary cells on processor
    vf_int bound_cells;

  public:
    // constructor
    General_Topology(const vf_int &, const vf_int &, const vf_int &, 
		     const std_string &);  
    
    // SERVICES (virtual from Topology)

    //! Return the number of global cells.
    int num_cells() const { return procs_per_cell.size(); }

    //! Return the number of cells on processor.
    int num_cells(int proc) const { return cells_per_proc[proc].size(); }

    //! Return the number of processors storing a global_cell.
    int num_procs(int gcell) const { return procs_per_cell[gcell-1].size(); } 

    // Get the global cell indices.
    inline int global_cell(int) const;
    inline int global_cell(int, int) const;

    // Get the local cell index.
    inline int local_cell(int) const;
    inline int local_cell(int, int) const;

    // Access boundary cell indices.
    inline int global_to_boundary(int, int) const;
    inline int boundary_to_global(int, int) const;

    //! Get the number of boundary cells on a given processor.
    inline int get_boundary_cells(int) const;  

    // Get a list of global cells on processor.
    inline sf_int get_cells(int) const;

    // Get a list of processors that a global cell is on.
    inline sf_int get_procs(int) const;

    // Diagnostics.
    void print(std::ostream &) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

//! Overloaded stream output operator.
std::ostream& operator<<(std::ostream &out, const Topology &object);

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS FOR GENERAL_TOPOLOGY
//---------------------------------------------------------------------------//
/*!
 * \brief Return the global cell index given a local cell index.
 *
 * This function is used on a processor to determine the global cell index
 * for a local cell index.
 *
 * The local cell index entered must be in the range of local cells that live
 * on the processor.  This range can be queried by calls to
 * Topology::num_cells with the current processor id as the argument.
 *
 * \param local_cell local cell index on a processor
 * \return global cell index
 */
int General_Topology::global_cell(int local_cell) const
{
    Require (local_cell > 0 &&
	     local_cell <= cells_per_proc[C4::node()].size());
    return cells_per_proc[C4::node()][local_cell-1];
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Return the global cell index given a local cell index and
 * processor id.
 *
 * This function is used on any processor to get a global cell index given a
 * local cell index and the processor id on which that local cell resides.
 * This function is identical to Topology::global_cell with the exception
 * that any processor can be entered.
 *
 * The local cell index entered must be in the range of local cells that live
 * on the processor id.  This range can be queried by calls to
 * Topology::num_cells with the processor id as the argument.
 *
 * \param local_cell local cell index on processor proc
 * \param proc processor on which local_cell resides
 * \return global cell index
 */
int General_Topology::global_cell(int local_cell, int proc) const
{
    Require (proc < C4::nodes());
    Require (local_cell > 0 && local_cell <= cells_per_proc[proc].size());
    return cells_per_proc[proc][local_cell-1];
}

//---------------------------------------------------------------------------//
/*!  
 * \brief Return the local cell index on the current processor given a global
 * cell index.
 *
 * This function is used on a processor to get a local cell index given a
 * global cell index.  If the requested global cell does not live on the
 * current processor a value of zero is returned.  Otherwise, the local cell
 * index for the particular global cell on the current processor is returned.
 *
 * \param global_cell requested global cell index
 * \return local cell index of global cell on processor or 0 if the global
 * cell does not live on the current processor
 */
int General_Topology::local_cell(int global_cell) const
{
    Require (global_cell > 0 && global_cell <= num_cells());

    // get the iterator location of the desired cell, we need to use const
    // iterators because that is what find returns-->see KAI reponse on
    // 19-JUN-98 for details
    int proc = C4::node();
    sf_int::const_iterator itr = std::find(cells_per_proc[proc].begin(), 
					   cells_per_proc[proc].end(), 
					   global_cell);

    // if the cell is here return it, remember the local_cell is [1,N]
    // whereas the dimensionality is [0,N-1]; thus we need to add one to the
    // iterator position to return the actual local_cell index
    if (itr != cells_per_proc[proc].end())
	return (itr - cells_per_proc[proc].begin()) + 1;

    // if the global_cell does not live on this node return 0
    return 0;
}


//---------------------------------------------------------------------------//
/*!  
 * \brief Return the local cell index given a global cell index and
 * processor id.
 *
 * This function is used on any processor to get a local cell index given a
 * global cell index and a processor id.  If the requested global cell does
 * not live on proc a value of zero is returned.  Otherwise, the local cell
 * index for the particular global cell on processor proc is returned.  It is
 * identical to Topology::local_cell with the exception that any processor
 * can be entered.
 *
 * \param global_cell requested global cell index
 * \param proc processor
 * \return local cell index of global cell on proc or 0 if the global cell
 * does not live on proc
 */
int General_Topology::local_cell(int global_cell, int proc) const
{
    Require (global_cell > 0 && global_cell <= num_cells());
    Require (proc < C4::nodes());

    // get the iterator location of the desired cell, we need to use const
    // iterators because that is what find returns-->see KAI reponse on
    // 19-JUN-98 for details
    sf_int::const_iterator itr = std::find(cells_per_proc[proc].begin(),
					   cells_per_proc[proc].end(), 
					   global_cell);

    // if the cell is here return it, remember the local_cell is [1,N]
    // whereas the dimensionality is [0,N-1]; thus we need to add one to the
    // iterator position to return the actual local_cell index
    if (itr != cells_per_proc[proc].end())
	return (itr - cells_per_proc[proc].begin()) + 1;

    // if the global_cell does not live on this node return 0
    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the boundary cell index for a global cell index and
 * processor id.
 *
 * This function gives the index of a boundary cell on processor for a given
 * global cell index.
 *
 * \param global_cell requested global cell index
 * \param proc processor id
 * \return boundary cell index of global cell on proc or 0 if the global cell
 * is not a boundary cell on the processor
 */
int General_Topology::global_to_boundary(int global_cell, int proc) const
{
    Require (global_cell > 0 && global_cell <= num_cells());
    Require (proc < C4::nodes());

    // leave if we have no boundary cells
    if (bound_cells.size() == 0) 
	return 0;
    else if (bound_cells[proc].size() == 0)
	return 0;

    // get the iterator location of the desired cell, we need to use const
    // iterators because that is what find returns-->see KAI reponse on
    // 19-JUN-98 for details
    sf_int::const_iterator itr = std::find(bound_cells[proc].begin(),
					   bound_cells[proc].end(), 
					   global_cell);

    // if the cell is here return it, remember the boundary cell is [1,N]
    // whereas the dimensionality is [0,N-1]; thus we need to add one to the
    // iterator position to return the actual boundary cell index
    if (itr != bound_cells[proc].end())
	return (itr - bound_cells[proc].begin()) + 1;

    // if the global_cell does not live on this node return 0
    return 0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the global cell index of a given boundary cell on processor.
 *
 * This function gives the index of a global cell for a given boundary cell
 * on processor.
 *
 * The boundary cell must be in the range of boundary cells on a processor.
 * This range may be queried by calls to get_boundary_cells() on a processor.
 *
 * \param boundary_cell boundary cell index, must be in range of boundary
 * cells on processor
 * \param proc processor id
 * \return global cell index for given boundary cell 
 */
int General_Topology::boundary_to_global(int boundary_cell, int proc) const
{
    Require (boundary_cell > 0 && boundary_cell <= bound_cells[proc].size()); 
    Require (proc < C4::nodes());

    return bound_cells[proc][boundary_cell-1];
}

//---------------------------------------------------------------------------//

int General_Topology::get_boundary_cells(int proc) const
{
    // if bound_cells.size is zero than this is a full replication topology
    if (bound_cells.size() == 0)
    {
	Check (get_parallel_scheme() == "replication");
	return 0;
    }
    
    return bound_cells[proc].size();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a local-global cell list on processor.
 *
 * This function returns a vector<int> of global cells on processor.  The
 * vector is dimensioned [0,N-1] where N is the number of local cells on the
 * processor.  This vector can be used to map local cells to global cells on
 * a processor.  The local cell index for each global cell is the vector
 * index + 1
 *
 * \param proc processor id
 * \return vector<int> of global cells on proc */
Topology::sf_int General_Topology::get_cells(int proc) const
{
    Require (proc < C4::nodes());
    return cells_per_proc[proc];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a list of processor that a global_cell lives on.
 *
 * This function returns a list of processors that have a copy of a queried
 * global cell.  It returns a vector<int> of processor ids, the vector is
 * dimensioned [0,N] where N+1 is the total number of processors that the
 * global cell lives on.
 *
 * \param global_cell queried global cell index
 * \return vector<int> of processor ids
 */
Topology::sf_int General_Topology::get_procs(int global_cell) const
{
    Require (global_cell > 0 && global_cell <= num_cells());
    return procs_per_cell[global_cell-1];
}

} // end namespace rtt_mc

#endif                          // __mc_General_Topology_hh__

//---------------------------------------------------------------------------//
//                              end of mc/General_Topology.hh
//---------------------------------------------------------------------------//
