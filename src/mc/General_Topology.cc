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

//===========================================================================//
// GENERAL_TOPOLOGY CLASS DEFINITIONS
//===========================================================================//

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
// PACK
//---------------------------------------------------------------------------//
/*!
 * \brief Pack up the General_Topology and return a smart pointer to a
 * Topology::Pack base class.

 * Data in the General_Topology is packed into a Topology::Pack base class
 * pointer.  The actual data resides in the General_Topology::Pack derived
 * class.  The pointer can be used to access pointer to the data for message
 * passing and output operations.

 */
Topology::SP_Pack General_Topology::pack() const
{
    // we are returning a smart pointer to the base class pack object
    Topology::SP_Pack return_pack;

    //>>> store the topology type

    // indicator for topology type
    int indicator = 0;
    if (get_parallel_scheme() == "replication")
	indicator = 1;
    else if (get_parallel_scheme() == "DD")
	indicator = 2;
    else if (get_parallel_scheme() == "DD/replication")
	indicator = 3;
    else
	Insist (0, "Invalid topology specified!");

    //>>> calculate the size of the data
    
    int size = 0;

    // size of cells_per_proc data
    for (int i = 0; i < cells_per_proc.size(); i++)
	size += cells_per_proc[i].size();

    // size of procs_per_cell data
    for (int i = 0; i < procs_per_cell.size(); i++)
	size += procs_per_cell[i].size();

    // size of boundary cell data
    for (int i = 0; i < bound_cells.size(); i++)
	size += bound_cells[i].size();

    // size indicators: size of each data element
    size += 3 + cells_per_proc.size() + procs_per_cell.size() +
	bound_cells.size();

    // allocate and integer array to hold everything
    int *data = new int[size];
    int ic    = 0;
    
    //>>> put sizes into array
    
    // cells_per_proc
    data[ic++] = cells_per_proc.size();
    for (int i = 0; i < cells_per_proc.size(); i++)
    {
	data[ic++] = cells_per_proc[i].size();
	for (int j = 0; j < cells_per_proc[i].size(); j++)
	    data[ic++] = cells_per_proc[i][j];
    }
    
    // procs_per_cell
    data[ic++] = procs_per_cell.size();
    for (int i = 0; i < procs_per_cell.size(); i++)
    {
	data[ic++] = procs_per_cell[i].size();
	for (int j = 0; j < procs_per_cell[i].size(); j++)
	    data[ic++] = procs_per_cell[i][j];
    }
    
    // bound_cells
    data[ic++] = bound_cells.size();
    for (int i = 0; i < bound_cells.size(); i++)
    {
	data[ic++] = bound_cells[i].size();
	for (int j = 0; j < bound_cells[i].size(); j++)
	    data[ic++] = bound_cells[i][j];
    }
    Check (ic == size);

    // make the pack object, no need to clean up data, the pack object will
    // do that
    return_pack = new General_Topology::Pack(indicator, size, data);
    return return_pack;
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

//===========================================================================//
// GENERAL_TOPOLOGY::PACK CLASS DEFINITION
//===========================================================================//
/*!
 * \brief Constructor.

 * Construct a General_Topology::Pack instance.  Once allocated integer data
 * is given to the General_Topology::Pack constructor in the form of an int*,
 * the Pack object owns it.  When the Pack object goes out of scope it will
 * delete clean up the memory.  In general, Pack objects are only created by
 * calling the General_Topology::pack() function.

 * \param i indicator of parallel topology: 1="replication"; 2="DD";
 * 3="DD/replication"
 * \param s size of integer data stream
 * \param d pointer to integer data stream

 */
General_Topology::Pack::Pack(int i, int s, int *d)
    : Topology::Pack(),
		indicator(i),
		data(d),
		size(s)
{
    // nothing more to do here
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.

 * Do copy construction while preserving memory.  This is not a reference
 * counted class so data is copied from one class to the other during
 * function calls and the like (wherever a copy constructor is called).

 */
General_Topology::Pack::Pack(const Pack &rhs)
    : Topology::Pack(),
		indicator(rhs.indicator),
		data(new int[rhs.size]),
		size(rhs.size)
{
    Require (indicator > 0 && indicator <= 3);

    // fill up new data array
    for (int i = 0; i < size; i++)
	data[i] = rhs.data[i];
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.

 * Cleans up memory when the Pack object goes out of scope.  Once allocated
 * pointers are given to the Pack object the Pack object takes control of
 * them.

 */
General_Topology::Pack::~Pack()
{
    delete [] data;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return parallel scheme description.

 * This is the descriptive string that matches the parallel scheme indicator: 
 * \arg 1 = "replication"
 * \arg 2 = "DD"
 * \arg 3 = "DD/replication"

 */
Topology::std_string General_Topology::Pack::get_parallel_scheme() const
{
    std_string desc;
    if (indicator == 1)
	desc = "replication";
    else if (indicator == 2)
	desc = "DD";
    else if (indicator == 3)
	desc == "DD/replication";
    else
	Insist (0, "Invalid indicator for the parallel scheme!");

    return desc;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpack a General_Topology::Pack object into a smart pointer to a
 * rtt_mc::Topology.

 * The unpack function creates a smart pointer (rtt_dsxx::SP) to the Topology
 * base class from the data stored inside the General_Topology::Pack object.
 * Obviously, the General_Topology::Pack class creates a Topology base class
 * pointer to a General_Topology object.

 * \return rtt_dsxx::SP to a Topology base class pointer that points to a
 * General_Topology derived class instance

 */
Topology::SP_Topology General_Topology::Pack::unpack() const
{
    // return a SP to a Topolgoy base class
    SP_Topology topology;

    // unpack the topology

    //>>> Get the topology type
    std_string desc = get_parallel_scheme();

    // counter
    int ic = 0;

    //>>> Build the requisite fields
    vf_int cpp;
    vf_int ppc;
    vf_int bc;

    // cells_per_proc
    cpp.resize(data[ic++]);
    for (int i = 0; i < cpp.size(); i++)
    {
	cpp[i].resize(data[ic++]);
	for (int j = 0; j < cpp[i].size(); j++)
	    cpp[i][j] = data[ic++];
    }

    // procs_per_cell
    ppc.resize(data[ic++]);
    for (int i = 0; i < ppc.size(); i++)
    {
	ppc[i].resize(data[ic++]);
	for (int j = 0; j < ppc[i].size(); j++)
	    ppc[i][j] = data[ic++];
    }

    // boundary cells
    bc.resize(data[ic++]);
    for (int i = 0; i < bc.size(); i++)
    {
	bc[i].resize(data[ic++]);
	for (int j = 0; j < bc[i].size(); j++)
	    bc[i][j] = data[ic++];
    }
    Check (ic == size);
    
    // rebuild the topology
    topology = new General_Topology(cpp, ppc, bc, desc);
    Ensure (topology);

    return topology;
}

} // end of rtt_mc

//---------------------------------------------------------------------------//
//                              end of General_Topology.cc
//---------------------------------------------------------------------------//
