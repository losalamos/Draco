//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/General_Topology.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 30 13:19:05 1999
 * \brief  General_Topology implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "General_Topology.hh"
#include <iomanip>

namespace rtt_mc
{

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::setw;
using std::setiosflags;

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
/*!
 * \brief General_Topology constructor.
 *
 * The General_Topology constructor builds a General_Topology instance given
 * complete state data about the topology.  We assume that some sort of
 * Topology builder class will be used to sift information from the mesh,
 * sources, and problem setup to appropriately build some sort of spatial
 * decomposition.
 *
 * \param cpp Topology::vf_int of global cells per processor.  The size of
 * cpp is the total number of processors in the problem.  The size of each
 * cpp[i] is the number of cells on processor i.  The return value of each
 * [processor][local_cell-1] query is the global cell index.
 *
 * \param ppc Topology::vf_int of processors per global cell.  The size of
 * ppc is the total number of cells in the problem.  The size of each ppc[i]
 * is the number of processors that store global cell i.  The return value of
 * each [global_cell-1][j] is a processor id that the global_cell lives on.
 *
 * \param bc Topology::vf_int of boundary cells on each processor.  The size
 * of bc is the number of processors in the problem. If ps != "replication".
 * If ps == "replication" than this entry must have size 0.  If ps !=
 * "replication" than the return value of each bc[processor][bound_cell-1] is
 * the global cell index.  The size of each bc[processor] is the number of
 * boundary cells on processor.
 *
 * \param ps string that describes the spatial partitioning.  Current options
 * are "replication", "DD", and "DD/replication".  
 */
General_Topology::General_Topology(const vf_int &cpp, const vf_int &ppc, 
				   const vf_int &bc, const std_string &ps)
    : Topology(ps), cells_per_proc(cpp), procs_per_cell(ppc), bound_cells(bc)
{
    int nodes = C4::nodes();
    Ensure (cells_per_proc.size() == nodes);

    // check (if we are at the appropriate DBC level) to make sure that
    // bound_cells.size() == 0 if we are doing full replication
#if DBC & 4
    if (get_parallel_scheme() == "replication")
    {
	Ensure (bound_cells.size() == 0);
    }
    else
    {
	Ensure (bound_cells.size() == nodes);
    }
#endif
}

//---------------------------------------------------------------------------//
// PRINT FOR DIAGNOSTICS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out data in the topology class for diagnostic purposes.
 */
void General_Topology::print(ostream &out) const
{
    out << endl;
    out << "*** GENERAL_TOPOLOGY DATA ***" << endl;
    out << "-----------------------------" << endl;

    out << endl;

    out << "Parallel Topology: " << get_parallel_scheme() << endl;
    out << "-----------------------------" << endl;

    // Print out cells_per_proc
    out << "Global Cells per Processor:" << endl;
    for (int i = 0; i < cells_per_proc.size(); i++)
    {
	out << setw(10) << setiosflags(ios::right) << "Node:" << setw(5) 
	    << i 
	    << setw(15) << setiosflags(ios::right) << "Local Cell" 
	    << setw(15) << setiosflags(ios::right) << "Global Cell"
	    << endl;
	for (int j = 0; j < cells_per_proc[i].size(); j++)
	    out << setw(30) << j+1 << setw(15) << cells_per_proc[i][j] 
		<< endl; 
    }

    out << endl;
}

} // end of rtt_mc

//---------------------------------------------------------------------------//
//                              end of General_Topology.cc
//---------------------------------------------------------------------------//
